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

#pragma once

#include <memory>
#include <map>
#include <Eigen/Core>

#include <RBDyn/FD.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/Joint.h>
#include <RBDyn/Coriolis.h>

#include <jrl-qp/experimental/BoxAndSingleConstraintSolver.h>

namespace torque_control
{

typedef std::map<std::string, int> ElapsedTimeMap;

class TorqueFeedbackTerm
{
public:
  enum TorqueControlType
  {
    None,
    IntegralTerm,
    PassivityPIDTerm
  };

  TorqueFeedbackTerm(const std::vector<rbd::MultiBody> & mbs,
                     int robotIndex,
                     const std::shared_ptr<rbd::ForwardDynamics> fd);

  const Eigen::VectorXd & P() const
  {
    return P_;
  }

  const Eigen::VectorXd& gammaD() const
  {
    return gammaD_;
  }

  virtual void computeTerm(const rbd::MultiBody & mb,
                           const rbd::MultiBodyConfig & mbc_real,
                           const rbd::MultiBodyConfig & mbc_calc) = 0;

  ElapsedTimeMap & getElapsedTimes();

 protected:

  int nrDof_;
  std::shared_ptr<rbd::ForwardDynamics> fd_;

  Eigen::VectorXd P_;
  Eigen::VectorXd gammaD_;

  Eigen::LLT<Eigen::MatrixXd> LLT_;

  ElapsedTimeMap elapsed_;

  void computeGammaD();
};

class IntegralTerm : public TorqueFeedbackTerm
{
 public:

  enum IntegralTermType
  {
    None,
    Simple,
    PassivityBased
  };

  enum VelocityGainType
  {
    Diagonal,
    MassDiagonal,
    MassMatrix
  };
  
  /// TorqueL is the lower bound for the torque
  /// TorqueU is the upper bound for the torque

  IntegralTerm(const std::vector<rbd::MultiBody> & mbs,
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
               double timeStep);

  void computeGain(const rbd::MultiBody & mb,
		   const rbd::MultiBodyConfig & mbc_real);

  void computeTerm(const rbd::MultiBody & mb,
                   const rbd::MultiBodyConfig & mbc_real,
                   const rbd::MultiBodyConfig & mbc_calc) override;

  void computeTerm(const rbd::MultiBody & mb,
                   const rbd::MultiBodyConfig & mbc_real,
                   const rbd::MultiBodyConfig & mbc_calc,
                   const Eigen::VectorXd & diff_torques);

  void enableNullSpaceCompliance(bool b);

  const Eigen::MatrixXd & CoriolisFactorization() const
  {
    return C_;
  }

  void setTaskSpaceJacobian(const std::string & descriptor, const Eigen::MatrixXd & Jac)
  {
    taskSpaceJacobians_[descriptor] = &Jac;
  }
  
  void setTaskSpaceJacobianDot(const std::string & descriptor, const Eigen::MatrixXd & JacDot)
  {
    taskSpaceJacobianDots_[descriptor] = &JacDot;
  }

  void setTaskSpaceErrorP(const std::string & descriptor, const Eigen::Vector3d & errorP)
  {
    taskSpaceErrorsP_[descriptor] = &errorP;
  }

  void setTaskSpaceErrorR(const std::string & descriptor, const Eigen::Matrix3d & errorR)
  {
    taskSpaceErrorsR_[descriptor] = &errorR;
  }
  
  void setTaskSpaceErrorV(const std::string & descriptor, const Eigen::Vector3d & errorV)
  {
    taskSpaceErrorsV_[descriptor] = &errorV;
  }

  void setTaskSpaceErrorW(const std::string & descriptor, const Eigen::Vector3d & errorW)
  {
    taskSpaceErrorsW_[descriptor] = &errorW;
  }
  
 protected:

  IntegralTermType intglTermType_;
  VelocityGainType velGainType_;
  double lambda_;

  bool enableNullSpaceCompliance_;

  rbd::Coriolis coriolis_;
  Eigen::MatrixXd C_;
  Eigen::VectorXd coriolisTerm_;
  Eigen::MatrixXd K_;

  Eigen::VectorXd intdjs_slow_;
  Eigen::VectorXd intdjs_fast_;
  Eigen::VectorXd previousS_;

  Eigen::VectorXd fastFilteredS_;
  Eigen::VectorXd slowFilteredS_;

  double phiSlow_;
  double phiFast_;
  double expPhiSlow_;
  double expPhiFast_;
  double fastFilterWeight_;

  Eigen::Vector3d maxLinAcc_, maxAngAcc_;
  Eigen::Vector6d curMaxFBWrench_, targetMaxFBWrench_;
  bool initializedMaxFBWrenches;

  Eigen::VectorXd torqueL_, torqueU_;
  double currentPerc_;
  double targetPerc_;

  int floatingBaseIndex_;
  jrl::qp::experimental::BoxAndSingleConstraintSolver solver_;

  double timeStep_;

  std::map<std::string, const Eigen::MatrixXd*> taskSpaceJacobians_;
  std::map<std::string, const Eigen::MatrixXd*> taskSpaceJacobianDots_;
  
  std::map<std::string, const Eigen::Vector3d*> taskSpaceErrorsP_;
  std::map<std::string, const Eigen::Matrix3d*> taskSpaceErrorsR_;
  std::map<std::string, const Eigen::Vector3d*> taskSpaceErrorsV_;
  std::map<std::string, const Eigen::Vector3d*> taskSpaceErrorsW_;
};

class PassivityPIDTerm : public TorqueFeedbackTerm
{
 public:

  PassivityPIDTerm(const std::vector<rbd::MultiBody> & mbs, int robotIndex,
                   const std::shared_ptr<rbd::ForwardDynamics> fd, double timeStep,
                   double beta, double lambda, double mu, double sigma, double cis);

  void computeTerm(const rbd::MultiBody & mb,
                   const rbd::MultiBodyConfig & mbc_real,
                   const rbd::MultiBodyConfig & mbc_calc) override;

  const Eigen::MatrixXd & CoriolisFactorization() const
  {
    return C_;
  }

 private:

  double dt_;
  double beta_, lambda_, mu_, sigma_, cis_;

  Eigen::VectorXd EPrev_;

  rbd::Coriolis coriolis_;
  Eigen::MatrixXd C_;

  Eigen::VectorXd errorParam(rbd::Joint::Type type,
                             const std::vector<double> & q_ref,
                             const std::vector<double> & q_hat);
};


} // namespace torque_control
