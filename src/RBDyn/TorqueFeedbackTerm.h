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
#include <Eigen/Core>

#include <RBDyn/FD.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <RBDyn/Joint.h>

namespace torque_control
{

class TorqueFeedbackTerm
{
 public:

  enum TorqueControlType
  {
    None,
    IntegralTerm,
    PassivityPIDTerm
  };
  
  TorqueFeedbackTerm(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                     const std::shared_ptr<rbd::ForwardDynamics> fd);
  
  const Eigen::VectorXd& P() const
  {
    return P_;
  }
  
  const Eigen::VectorXd& gammaD() const
  {
    return gammaD_;
  }

  virtual void computeTerm(const rbd::MultiBody& mb,
                           const rbd::MultiBodyConfig& mbc_real,
                           const rbd::MultiBodyConfig& mbc_calc) = 0;

 protected:

  int nrDof_;
  std::shared_ptr<rbd::ForwardDynamics> fd_;
  
  Eigen::VectorXd P_;
  Eigen::VectorXd gammaD_;
  
  Eigen::LLT<Eigen::MatrixXd> LLT_;

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

  IntegralTerm(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
               const std::shared_ptr<rbd::ForwardDynamics> fd,
               IntegralTermType intglTermType, VelocityGainType velGainType,
               double lambda);
  
  void computeTerm(const rbd::MultiBody& mb,
                   const rbd::MultiBodyConfig& mbc_real,
                   const rbd::MultiBodyConfig& mbc_calc) override;
    
 private:
    
  IntegralTermType intglTermType_;
  VelocityGainType velGainType_;
  double lambda_;
};


class PassivityPIDTerm : public TorqueFeedbackTerm
{
 public:

  PassivityPIDTerm(const std::vector<rbd::MultiBody>& mbs, int robotIndex,
                   const std::shared_ptr<rbd::ForwardDynamics> fd, double timeStep,
                   double beta, double lambda, double mu, double sigma, double cis);
  
  void computeTerm(const rbd::MultiBody& mb,
                   const rbd::MultiBodyConfig& mbc_real,
                   const rbd::MultiBodyConfig& mbc_calc) override;
  
 private:

  double dt_;
  double beta_, lambda_, mu_, sigma_, cis_;

  Eigen::VectorXd EPrev_;

  Eigen::VectorXd errorParam(rbd::Joint::Type type,
                             std::vector<double> q_ref,
                             std::vector<double> q_hat);
};

 
} // namespace torque_control