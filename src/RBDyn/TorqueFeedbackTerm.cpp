// Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
//
// This file is part of RBDyn.
//
// RBDyn is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RBDyn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with RBDyn.  If not, see <http://www.gnu.org/licenses/>.

#include "RBDyn/TorqueFeedbackTerm.h"

#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

#include <time.h>

namespace torque_control
{

/**
 *    TorqueFeedbackTerm
 */

TorqueFeedbackTerm::TorqueFeedbackTerm(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                                       const std::shared_ptr<rbd::ForwardDynamics> fd)
  : nrDof_(mbs[robotIndex].nrDof()),
    fd_(fd),
    P_(Eigen::VectorXd::Zero(nrDof_)),
    gammaD_(Eigen::VectorXd::Zero(nrDof_))
{
  elapsed_ = {{"computeFbTerm-Gain", 0},
              {"computeFbTerm-Gain-Coriolis", 0},
              {"computeFbTerm-GammaD", 0}};
}

void TorqueFeedbackTerm::computeGammaD()
{
  // Alternative method to do
  // gammaD_ = fd_->H().inverse() * P_;

  gammaD_ = P_;
  LLT_.compute(fd_->H());
  LLT_.matrixL().solveInPlace(gammaD_);
  LLT_.matrixL().transpose().solveInPlace(gammaD_);
}

ElapsedTimeMap & TorqueFeedbackTerm::getElapsedTimes()
{
  return elapsed_;
}

/**
 *    IntegralTerm
 */

IntegralTerm::IntegralTerm(const std::vector<rbd::MultiBody> & mbs,
                           int robotIndex,
                           const std::shared_ptr<rbd::ForwardDynamics> fd,
                           IntegralTermType intglTermType,
                           VelocityGainType velGainType,
                           double lambda,
                           double perc,
                           const Eigen::Vector3d & maxLinAcc,
                           const Eigen::Vector3d & maxAngAcc,
                           const Eigen::VectorXd & torqueL,
                           const Eigen::VectorXd & torqueU,
                           double phiSlow,
                           double phiFast,
                           double fastFilterWeight,
                           double timeStep)
: TorqueFeedbackTerm(mbs, robotIndex, fd), intglTermType_(intglTermType), velGainType_(velGainType), lambda_(lambda),
  coriolis_(mbs[robotIndex]), C_(Eigen::MatrixXd::Zero(nrDof_, nrDof_)), coriolisTerm_(Eigen::VectorXd::Zero(nrDof_)),
  K_(Eigen::MatrixXd::Zero(nrDof_, nrDof_)), previousS_(Eigen::VectorXd::Zero(nrDof_)),
  fastFilteredS_(Eigen::VectorXd::Zero(nrDof_)), slowFilteredS_(Eigen::VectorXd::Zero(nrDof_)), phiSlow_(phiSlow),
  phiFast_(phiFast), expPhiSlow_(exp(-timeStep * phiSlow)), expPhiFast_(exp(-timeStep * phiFast)),
  fastFilterWeight_(fastFilterWeight), maxLinAcc_(maxLinAcc), maxAngAcc_(maxAngAcc), 
  curMaxFBWrench_(Eigen::Vector6d::Zero()), targetMaxFBWrench_(Eigen::Vector6d::Zero()), 
  initializedMaxFBWrenches(false),  timeStep_(timeStep), targetPerc_(perc), currentPerc_(perc),
   torqueL_(torqueL), torqueU_(torqueU), floatingBaseIndex_(-1),  solver_(mbs[robotIndex].nrJoints())
{  
  for(int i = 0; i < mbs[robotIndex].nrJoints(); i++)
  {
    if(mbs[robotIndex].joint(i).type() == rbd::Joint::Free)
    {
      floatingBaseIndex_ = mbs[robotIndex].jointPosInDof(i);

      /// set the floating base limits to infinity by default
      /// they will be replaced with wrench limit defined by the max acceleraTion
      /// this is useful for the transition
      torqueL_.segment<6>(floatingBaseIndex_).setConstant(- std::numeric_limits<double>::infinity());
      torqueU_.segment<6>(floatingBaseIndex_).setConstant( std::numeric_limits<double>::infinity());

      break;
    }
  }
}

void IntegralTerm::computeTerm(const rbd::MultiBody & mb,
                               const rbd::MultiBodyConfig & mbc_real,
                               const rbd::MultiBodyConfig & mbc_calc)
{
  if(intglTermType_ == Simple || intglTermType_ == PassivityBased)
  {
    computeGain(mb, mbc_real);

    Eigen::VectorXd alphaVec_ref = rbd::dofToVector(mb, mbc_calc.alpha);
    Eigen::VectorXd alphaVec_hat = rbd::dofToVector(mb, mbc_real.alpha);

    Eigen::VectorXd newS = alphaVec_ref - alphaVec_hat;

    slowFilteredS_ = expPhiSlow_ * slowFilteredS_ + newS - previousS_;
    fastFilteredS_ = expPhiFast_ * fastFilteredS_ + newS - previousS_;

    previousS_ = newS;

    Eigen::VectorXd filteredS = fastFilterWeight_ * fastFilteredS_ + (1 - fastFilterWeight_) * slowFilteredS_;

    
    ///compute the max torque allowed for the integral term
    Eigen::VectorXd torqueU_prime = torqueU_ * currentPerc_;
    Eigen::VectorXd torqueL_prime = torqueL_ * currentPerc_;

    if(floatingBaseIndex_ != -1)
    {
      if (!initializedMaxFBWrenches)
      {
        Eigen::Vector6d acc;
        acc << maxAngAcc_, maxLinAcc_;

        targetMaxFBWrench_ = curMaxFBWrench_ =
            fd_->H().block<6, 6>(floatingBaseIndex_, floatingBaseIndex_).diagonal().array() * acc.array().abs();

        initializedMaxFBWrenches = true;
      }
      torqueU_prime.segment<6>(floatingBaseIndex_) = curMaxFBWrench_;
      torqueL_prime.segment<6>(floatingBaseIndex_) = -curMaxFBWrench_;
    }

    /// make current limit converge to tarrget value
    /// this exponential convergence allows to avoid discontinuous torques
    currentPerc_ = expPhiSlow_ * currentPerc_ + (1 - expPhiSlow_) * targetPerc_;
    curMaxFBWrench_ = expPhiSlow_ * curMaxFBWrench_ + (1 - expPhiSlow_) * targetMaxFBWrench_;

    if(intglTermType_ == PassivityBased)
    {
      coriolisTerm_ = C_ * filteredS;

      /// CANCELED /// substract the Coriolis term to make sure that the whole integral term
      /// CANCELED /// is respecting tha antiwindup
      // torqueU_prime -= coriolisTerm_;
      // torqueL_prime -= coriolisTerm_;
    }

    P_ = K_ * filteredS;

    /// get the multiplier of the bound violation
    double epsilonU = (P_.array() / torqueU_prime.array()).maxCoeff();
    double epsilonL = (P_.array() / torqueL_prime.array()).maxCoeff();
    double epsilon = std::max(std::max(epsilonU, epsilonL), 1.);

    Eigen::MatrixXd serror(3, mb.nrDof());
    serror.row(0) = filteredS.transpose();
    serror.row(1) = alphaVec_hat.transpose();
    serror.row(2) = alphaVec_ref.transpose();

    std::cout << "Mehdi  serror   " << std::endl;
    std::cout << serror << std::endl;

    if(epsilon > 1)
    {
      double dotprod = P_.dot(filteredS) / epsilon;
      std::cout << "Mehdi QP Started" << std::endl;

      if(solver_.solve(P_, filteredS, dotprod, torqueL_prime, torqueU_prime) == jrl::qp::TerminationStatus::SUCCESS)
      {
        std::cout << "Mehdi Succeeded" << (solver_.solution() - P_).norm() << " " << (P_ / epsilon - P_).norm()
                  << std::endl;
        P_ = solver_.solution();
      }
      else
      {
        std::cout << "Mehdi QP FAILED" << std::endl;
        std::cerr << "Mehdi QP FAILED" << std::endl;

        P_ /= epsilon;
      }

      /// -----------------------------------------testing -------------------------------------
      {
        double epsilonU = (P_.array() / torqueU_prime.array()).maxCoeff();
        double epsilonL = (P_.array() / torqueL_prime.array()).maxCoeff();
        double epsilon = std::max(std::max(epsilonU, epsilonL), 1.);

        if(epsilon > 1 + 1e-5)
        {
          std::cout << "Mehdi Bounds not respected !" << std::endl;
          std::cerr << "Mehdi Bounds not respected !" << std::endl;
        }
        else
        {
          if(P_.dot(filteredS) < dotprod + 1e-5)
          {
            std::cout << "Mehdi dotProd not respected ! ref " << dotprod << " solution " << P_.dot(filteredS) << std::endl;
            std::cerr << "Mehdi dotProd not respected ! ref " << dotprod << " solution " << P_.dot(filteredS) << std::endl;
          }
        }
      }
    }

    if(intglTermType_ == PassivityBased)
    {
      P_ += coriolisTerm_;
    }

    computeGammaD();
  }
}

void IntegralTerm::computeTerm(const rbd::MultiBody & mb,
                               const rbd::MultiBodyConfig & mbc_real,
                               const rbd::MultiBodyConfig & mbc_calc,
                               const Eigen::VectorXd & diff_torques)
{
  computeTerm(mb, mbc_real, mbc_calc);

  std::cout << "Rafa, in IntegralTerm::computeTerm for smooth transition, diff_torques = " << diff_torques.transpose()
            << std::endl;

  
  if(intglTermType_ == Simple || intglTermType_ == PassivityBased)
  {
    Eigen::VectorXd diff_torques_for_antiwindup = diff_torques - coriolisTerm_;
    double percU = (diff_torques_for_antiwindup.array() / torqueU_.array()).maxCoeff();
    double percL = (diff_torques_for_antiwindup.array() / torqueL_.array()).maxCoeff();
    currentPerc_ = std::max(std::max(percU, percL), currentPerc_);
    if(floatingBaseIndex_!=-1)
    {
      curMaxFBWrench_ = curMaxFBWrench_.cwiseMax(diff_torques_for_antiwindup.segment<6>(floatingBaseIndex_).cwiseAbs());
    }      

    Eigen::MatrixXd L;   
    if(intglTermType_ == Simple)
    {
      L = K_;
    }
    else
    {
      L = K_ + C_;
    }

    if(fastFilterWeight_ < 1)
    {
      slowFilteredS_ = L.inverse() * diff_torques / (1 - fastFilterWeight_);
      fastFilteredS_.setZero();
    }
    else
    {
      slowFilteredS_.setZero();
      fastFilteredS_ = L.inverse() * diff_torques / fastFilterWeight_;
    }

     P_ = diff_torques;

    std::cout << "Rafa, in IntegralTerm::computeTerm for smooth transition, P_ = " << P_.transpose() << std::endl;
  }
}

void IntegralTerm::computeGain(const rbd::MultiBody & mb, const rbd::MultiBodyConfig & mbc_real)
{
  // std::cout << "Rafa, in IntegralTerm::computeTerm, fd_->H().rows() = " << fd_->H().rows() << std::endl;

  if(velGainType_ == MassMatrix)
  {
    K_ = lambda_ * fd_->H();
  }
  else if(velGainType_ == MassDiagonal)
  {
    K_ = lambda_ * fd_->H().diagonal().asDiagonal();
  }
  else
  {
    K_ = lambda_ * Eigen::MatrixXd::Identity(nrDof_, nrDof_);
  }

  clock_t time;

  if(intglTermType_ == PassivityBased)
  {
    time = clock();
    C_ = coriolis_.coriolis(mb, mbc_real);
    elapsed_.at("computeFbTerm-Gain-Coriolis") = (int)(clock() - time);
  }
}

/**
  *    PassivityPIDTerm
  */

PassivityPIDTerm::PassivityPIDTerm(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
                                   const std::shared_ptr<rbd::ForwardDynamics> fd, double timeStep,
                                   double beta, double lambda, double mu, double sigma, double cis)
  : TorqueFeedbackTerm(mbs, robotIndex, fd),
    dt_(timeStep),
    beta_(beta),
    lambda_(lambda),
    mu_(mu),
    sigma_(sigma),
    cis_(cis),
    EPrev_(Eigen::VectorXd::Zero(nrDof_)),
    coriolis_(mbs[robotIndex]),
    C_(Eigen::MatrixXd::Zero(nrDof_, nrDof_))
{
}

void PassivityPIDTerm::computeTerm(const rbd::MultiBody & mb,
                                   const rbd::MultiBodyConfig & mbc_real,
                                   const rbd::MultiBodyConfig & mbc_calc)
{
  const Eigen::MatrixXd & M = fd_->H();

  C_ = coriolis_.coriolis(mb, mbc_real);

  Eigen::MatrixXd Ka = beta_  * M;
  //Eigen::MatrixXd L  = sigma_ * M.diagonal().asDiagonal();
  Eigen::MatrixXd L  = sigma_ * M;

  Eigen::MatrixXd Kv = lambda_ * M + C_ + Ka;
  Eigen::MatrixXd Kp = mu_ * M + lambda_ * (C_ + Ka) + L;
  Eigen::MatrixXd Ki = mu_ * (C_ + Ka) + cis_ * lambda_ * L;

  Eigen::VectorXd alphaVec_ref = rbd::dofToVector(mb, mbc_calc.alpha);
  Eigen::VectorXd alphaVec_hat = rbd::dofToVector(mb, mbc_real.alpha);

  Eigen::VectorXd s = alphaVec_ref - alphaVec_hat;

  const std::vector<rbd::Joint>& joints = mb.joints();

  Eigen::VectorXd e = Eigen::VectorXd::Zero(nrDof_);

  int pos = 0;
  for (std::size_t i = 0; i < joints.size(); i++) {
    Eigen::VectorXd ei = errorParam(joints[i].type(), mbc_calc.q[i], mbc_real.q[i]);
    e.segment(pos, ei.size()) = ei;
    pos += ei.size();
  }

  Eigen::VectorXd E = EPrev_ + e * dt_;
  EPrev_ = E;

  P_ = Kv * s + Kp * e + Ki * E;

  computeGammaD();
}

Eigen::VectorXd PassivityPIDTerm::errorParam(rbd::Joint::Type type,
                                             const std::vector<double> & q_ref,
                                             const std::vector<double> & q_hat)
{
  Eigen::VectorXd e;

  switch (type) {

  case rbd::Joint::Rev:
  case rbd::Joint::Prism:

    e.resize(1);
    e[0] = (q_ref[0] - q_hat[0]);
    break;

  case rbd::Joint::Free:

    e.resize(6);

    Eigen::Quaterniond quat_ref(q_ref[0], q_ref[1], q_ref[2], q_ref[3]);
    Eigen::Quaterniond quat_hat(q_hat[0], q_hat[1], q_hat[2], q_hat[3]);

    Eigen::MatrixXd Re = quat_ref.normalized().toRotationMatrix() * quat_hat.normalized().toRotationMatrix().transpose();
    Eigen::MatrixXd Sw = Re.log();

    e[0] = Sw(2, 1);
    e[1] = Sw(0, 2);
    e[2] = Sw(1, 0);

    for (std::size_t i = 0; i < 3; i++)
      e[3 + i] = q_ref[4 + i] - q_hat[4 + i];

    break;
  }

  return e;
}

} // namespace torque_control
