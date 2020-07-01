/*
 * Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// associated header
#include "RBDyn/FD.h"

// includes
// RBDyn
#include "RBDyn/MultiBody.h"
#include "RBDyn/MultiBodyConfig.h"

// #include <iostream>

namespace rbd
{

ForwardDynamics::ForwardDynamics(const MultiBody & mb)
: H_(mb.nrDof(), mb.nrDof()), C_(mb.nrDof()), I_st_(mb.nrBodies()), F_(mb.nrJoints()), HIr_(mb.nrDof(), mb.nrDof()), 
  acc_(mb.nrBodies()), f_(mb.nrBodies()), tmpFd_(mb.nrDof()), dofPos_(mb.nrJoints()), ldlt_(mb.nrDof())
{
  HIr_.setZero();
  int dofP = 0;
  for(int i = 0; i < mb.nrJoints(); ++i)
  {
    F_[i].resize(6, mb.joint(i).dof());
    dofPos_[i] = dofP;
    dofP += mb.joint(i).dof();
  }
}

void ForwardDynamics::forwardDynamics(const MultiBody & mb, MultiBodyConfig & mbc)
{
  computeH(mb, mbc);
  computeC(mb, mbc);

  paramToVector(mbc.jointTorque, tmpFd_);
  ldlt_.compute(H_);
  tmpFd_ = ldlt_.solve(tmpFd_ - C_);

  vectorToParam(tmpFd_, mbc.alphaD);
}

void ForwardDynamics::computeHIr(const MultiBody& mb)
{
  for(int i = 0; i < mb.nrJoints(); ++i)
  {
    if(mb.joint(i).type() == Joint::Rev)
    {
      double gr = mb.joint(i).gearRatio();
      HIr_(dofPos_[i], dofPos_[i]) = mb.joint(i).rotorInertia() * gr * gr;
    }
  }
}

void ForwardDynamics::computeH(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  const std::vector<Body> & bodies = mb.bodies();
  const std::vector<Joint> & joints = mb.joints();
  const std::vector<int> & pred = mb.predecessors();

  H_.setZero();
  for(std::size_t i = 0; i < bodies.size(); ++i)
  {
    I_st_[i] = bodies[i].inertia();
  }

  for(int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i)
  {
    if(pred[i] != -1)
    {
      const sva::PTransformd & X_p_i = mbc.parentToSon[i];
      I_st_[pred[i]] += X_p_i.transMul(I_st_[i]);
    }

    for(int dof = 0; dof < joints[i].dof(); ++dof)
    {
      F_[i].col(dof).noalias() = (I_st_[i] * sva::MotionVecd(mbc.motionSubspace[i].col(dof))).vector();
    }

    if ( joints[i].dof() > 0)
    {
      H_.block(dofPos_[i], dofPos_[i], joints[i].dof(), joints[i].dof()).noalias() =
        mbc.motionSubspace[i].transpose() * F_[i];
    }

    int j = i;
    while(pred[j] != -1)
    {
      const sva::PTransformd & X_p_j = mbc.parentToSon[j];
      for(int dof = 0; dof < joints[i].dof(); ++dof)
      {
        F_[i].col(dof) = X_p_j.transMul(sva::ForceVecd(F_[i].col(dof))).vector();
      }
      j = pred[j];

      if(joints[j].dof() != 0)
      {
        H_.block(dofPos_[i], dofPos_[j], joints[i].dof(), joints[j].dof()).noalias() =
            F_[i].transpose() * mbc.motionSubspace[j];

        H_.block(dofPos_[j], dofPos_[i], joints[j].dof(), joints[i].dof()).noalias() =
            H_.block(dofPos_[i], dofPos_[j], joints[i].dof(), joints[j].dof()).transpose();
      }
    }
  }

  H_.noalias() = H_ + HIr_;
}

void ForwardDynamics::computeC(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  const std::vector<Body> & bodies = mb.bodies();
  const std::vector<Joint> & joints = mb.joints();
  const std::vector<int> & pred = mb.predecessors();

  sva::MotionVecd a_0(Eigen::Vector3d::Zero(), mbc.gravity);

  for(std::size_t i = 0; i < bodies.size(); ++i)
  {
    const sva::PTransformd & X_p_i = mbc.parentToSon[i];

    const sva::MotionVecd & vj_i = mbc.jointVelocity[i];

    const sva::MotionVecd & vb_i = mbc.bodyVelB[i];

    if(pred[i] != -1)
      acc_[i] = X_p_i * acc_[pred[i]] + vb_i.cross(vj_i);
    else
      acc_[i] = X_p_i * a_0 + vb_i.cross(vj_i);

    f_[i] = bodies[i].inertia() * acc_[i] + vb_i.crossDual(bodies[i].inertia() * vb_i)
            - mbc.bodyPosW[i].dualMul(mbc.force[i]);
  }

  for(int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i)
  {
    if (joints[i].dof() > 0)
    {
      C_.segment(dofPos_[i], joints[i].dof()).noalias() = mbc.motionSubspace[i].transpose() * f_[i].vector();
    }    

    if(pred[i] != -1)
    {
      const sva::PTransformd & X_p_i = mbc.parentToSon[i];
      f_[pred[i]] += X_p_i.transMul(f_[i]);
    }
  }
}

void ForwardDynamics::sForwardDynamics(const MultiBody & mb, MultiBodyConfig & mbc)
{
  checkMatchParentToSon(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);
  checkMatchJointVelocity(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchBodyPos(mb, mbc);
  checkMatchForce(mb, mbc);
  checkMatchJointTorque(mb, mbc);

  checkMatchAlphaD(mb, mbc);

  forwardDynamics(mb, mbc);
}

void ForwardDynamics::sComputeH(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchParentToSon(mb, mbc);
  // checkMatchMotionSubspace(mb, mbc);  // Commented out by Rafa as a test

  computeH(mb, mbc);
}

void ForwardDynamics::sComputeC(const MultiBody & mb, const MultiBodyConfig & mbc)
{
  checkMatchParentToSon(mb, mbc);
  checkMatchMotionSubspace(mb, mbc);
  checkMatchJointVelocity(mb, mbc);
  checkMatchBodyVel(mb, mbc);
  checkMatchBodyPos(mb, mbc);
  checkMatchForce(mb, mbc);

  computeC(mb, mbc);
}

/*
Eigen::Matrix3d ForwardDynamics::SkewSymmetric(const Eigen::Vector3d& v)
{
  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  R(0,1) = -(R(1,0) = v[2]);
  R(2,0) = -(R(0,2) = v[1]);
  R(1,2) = -(R(2,1) = v[0]);
  
  return R;
}
*/

} // namespace rbd
