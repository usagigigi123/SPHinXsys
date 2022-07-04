/**
 * @file 	fluid_dynamics.cpp
 * @author	Chi ZHang and Xiangyu Hu
 */

#include "fluid_dynamics_inner.h"
#include "fluid_dynamics_inner.hpp"

namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		FluidInitialCondition::
			FluidInitialCondition(FluidBody &fluid_body)
			: ParticleDynamicsSimple(fluid_body), FluidDataSimple(fluid_body),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_) {}
		//=================================================================================================//
		DensitySummationInner::DensitySummationInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamicsWithUpdate(*inner_relation.sph_body_),
			  FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), mass_(particles_->mass_),
			  rho_sum_(particles_->rho_sum_),
			  W0_(sph_adaptation_->getKernel()->W0(Vecd(0))),
			  rho0_(particles_->rho0_), inv_sigma0_(1.0 / particles_->sigma0_) {}
		//=================================================================================================//
		void DensitySummationInner::Interaction(size_t index_i, Real dt)
		{
			/** Inner interaction. */
			Real sigma = W0_;
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
		}
		//=================================================================================================//
		void DensitySummationInner::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] = ReinitializedDensity(rho_sum_[index_i], rho0_, rho_n_[index_i]);
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
		}
				//=================================================================================================//
		UpdateViscosity::UpdateViscosity(BaseBodyRelationInner &inner_relation, 
		            Real mu_0, Real lambda, Real a, Real b)
			: InteractionDynamics(*inner_relation.sph_body_),
			  FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), vel_n_(particles_->vel_n_),
			  mu_(material_->ReferenceViscosity()),
			  mu_0_(mu_0), lambda_(lambda), a_(a), b_(b),
			  shear_rate_(DynamicCast<ShearThinningFluidParticles>(this, sph_body_->base_particles_)->shear_rate_),
			  mu_shear_(DynamicCast<ShearThinningFluidParticles>(this, sph_body_->base_particles_)->mu_shear_)
		{}
		//=================================================================================================//
		void UpdateViscosity::Interaction(size_t index_i, Real dt)
		{
			const Vecd vel_i = vel_n_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];

			Matd shear_rate_tensor(0);
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				//viscous force
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
				Matd velocity_gradient_i = -SimTK::outer((vel_i - vel_n_[index_j]), nablaW_ij) * Vol_[index_j];
				shear_rate_tensor += 0.5 * (velocity_gradient_i + ~velocity_gradient_i);
			}
			Real shear_rate = sqrt(2.0) * shear_rate_tensor.norm();
			shear_rate_[index_i] = shear_rate;
			//mu_shear_[index_i] = mu_ + (mu_0_ - mu_) * (1 + log(1 + lambda_ * shear_rate)) / (1 + lambda_ * shear_rate);
			mu_shear_[index_i] = mu_ + (mu_0_ - mu_) / std::pow(1 + std::pow(lambda_ * shear_rate_[index_i], b_), a_);
		}
		//=================================================================================================//
		ViscousAccelerationInner::ViscousAccelerationInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_),
			  vel_n_(particles_->vel_n_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_),
			  mu_(material_->ReferenceViscosity()),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		void ViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			const Vecd &vel_i = vel_n_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];

				//viscous force
				vel_derivative = (vel_i - vel_n_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}

			dvel_dt_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		void ShearThinningViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			const Vecd &vel_i = vel_n_[index_i];
			Real mu_i = mu_shear_[index_i];

			Vecd acceleration(0), vel_derivative(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				//viscous force
				Real mu_j = mu_shear_[index_j];
				Real mu_ij = 2 * mu_i * mu_j / (mu_i + mu_j);
				vel_derivative = (vel_i - vel_n_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
				acceleration += 2.0 * mu_ij * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
				//acceleration += 2.0 * mu_i * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
			}

			dvel_dt_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		void AngularConservativeViscousAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];
			const Vecd &vel_i = vel_n_[index_i];

			Vecd acceleration(0);
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd &e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];

				/** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that 
				 * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
				Real v_r_ij = dot(vel_i - vel_n_[index_j], r_ij * e_ij);
				Real eta_ij = 8.0 * mu_ * v_r_ij / (r_ij * r_ij + 0.01 * smoothing_length_);
				acceleration += eta_ij * Vol_[index_j] / rho_i * inner_neighborhood.dW_ij_[n] * e_ij;
			}

			dvel_dt_prior_[index_i] += acceleration;
		}
		//=================================================================================================//
		TransportVelocityCorrectionInner::
			TransportVelocityCorrectionInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_),
			  pos_n_(particles_->pos_n_),
			  surface_indicator_(particles_->surface_indicator_), p_background_(0) {}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::setupDynamics(Real dt)
		{
			Real speed_max = particles_->speed_max_;
			Real density = material_->ReferenceDensity();
			p_background_ = 7.0 * density * speed_max * speed_max;
		}
		//=================================================================================================//
		void TransportVelocityCorrectionInner::Interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_n_[index_i];

			Vecd acceleration_trans(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				//acceleration for transport velocity
				acceleration_trans -= 2.0 * p_background_ * Vol_[index_j] * nablaW_ij / rho_i;
			}

			if (surface_indicator_[index_i] == 0)
				pos_n_[index_i] += acceleration_trans * dt * dt * 0.5;
		}
		//=================================================================================================//
		AcousticTimeStepSize::AcousticTimeStepSize(FluidBody &fluid_body)
			: ParticleDynamicsReduce<Real, ReduceMax>(fluid_body),
			  FluidDataSimple(fluid_body), rho_n_(particles_->rho_n_),
			  p_(particles_->p_), vel_n_(particles_->vel_n_),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength())
		{
			initial_reference_ = 0.0;
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return material_->getSoundSpeed(p_[index_i], rho_n_[index_i]) + vel_n_[index_i].norm();
		}
		//=================================================================================================//
		Real AcousticTimeStepSize::OutputResult(Real reduced_value)
		{
			particles_->signal_speed_max_ = reduced_value;
			//since the particle does not change its configuration in pressure relaxation step
			//I chose a time-step size according to Eulerian method
			return 0.6 * smoothing_length_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSize::AdvectionTimeStepSize(FluidBody &fluid_body, Real U_max)
			: ParticleDynamicsReduce<Real, ReduceMax>(fluid_body),
			  FluidDataSimple(fluid_body), vel_n_(particles_->vel_n_),
			  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength())
		{
			Real rho_0 = material_->ReferenceDensity();
			Real mu = material_->ReferenceViscosity();
			Real viscous_speed = mu / rho_0 / smoothing_length_;
			Real u_max = SMAX(viscous_speed, U_max);
			initial_reference_ = u_max * u_max;
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::ReduceFunction(size_t index_i, Real dt)
		{
			return vel_n_[index_i].normSqr();
		}
		//=================================================================================================//
		Real AdvectionTimeStepSize::OutputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			particles_->speed_max_ = speed_max;
			return 0.25 * smoothing_length_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		AdvectionTimeStepSizeForImplicitViscosity::
			AdvectionTimeStepSizeForImplicitViscosity(FluidBody &fluid_body, Real U_max)
			: AdvectionTimeStepSize(fluid_body, U_max)
		{
			initial_reference_ = U_max * U_max;
		}
		//=================================================================================================//
		VorticityInner::
			VorticityInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), vel_n_(particles_->vel_n_)
		{
			particles_->registerAVariable(vorticity_, "VorticityInner");
			particles_->addAVariableToWrite<AngularVecd>("VorticityInner");
		}
		//=================================================================================================//
		void VorticityInner::Interaction(size_t index_i, Real dt)
		{
			const Vecd &vel_i = vel_n_[index_i];

			AngularVecd vorticity(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd r_ij = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

				Vecd vel_diff = vel_i - vel_n_[index_j];
				vorticity += SimTK::cross(vel_diff, r_ij) * Vol_[index_j] * inner_neighborhood.dW_ij_[n];
			}

			vorticity_[index_i] = vorticity;
		}
		//=================================================================================================//
		BaseRelaxation::BaseRelaxation(BaseBodyRelationInner &inner_relation)
			: ParticleDynamics1Level(*inner_relation.sph_body_),
			  FluidDataInner(inner_relation),
			  Vol_(particles_->Vol_), mass_(particles_->mass_), rho_n_(particles_->rho_n_),
			  p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
			  dvel_dt_(particles_->dvel_dt_),
			  dvel_dt_prior_(particles_->dvel_dt_prior_) {}
		//=================================================================================================//
		BasePressureRelaxation::
			BasePressureRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BasePressureRelaxation::Initialization(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
			Vol_[index_i] = mass_[index_i] / rho_n_[index_i];
			p_[index_i] = material_->getPressure(rho_n_[index_i]);
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void BasePressureRelaxation::Update(size_t index_i, Real dt)
		{
			vel_n_[index_i] += dvel_dt_[index_i] * dt;
		}
		//=================================================================================================//
		Vecd BasePressureRelaxation::computeNonConservativeAcceleration(size_t index_i)
		{
			Real rho_i = rho_n_[index_i];
			Real p_i = p_[index_i];
			Vecd acceleration = dvel_dt_prior_[index_i];
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real dW_ij = inner_neighborhood.dW_ij_[n];
				const Vecd &e_ij = inner_neighborhood.e_ij_[n];

				Real rho_j = rho_n_[index_j];
				Real p_j = p_[index_j];

				Real p_star = (rho_i * p_j + rho_j * p_i) / (rho_i + rho_j);
				acceleration += (p_i - p_star) * Vol_[index_j] * dW_ij * e_ij / rho_i;
			}
			return acceleration;
		}
		//=================================================================================================//
		BaseDensityRelaxation::
			BaseDensityRelaxation(BaseBodyRelationInner &inner_relation) : BaseRelaxation(inner_relation) {}
		//=================================================================================================//
		void BaseDensityRelaxation::Initialization(size_t index_i, Real dt)
		{
			pos_n_[index_i] += vel_n_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void BaseDensityRelaxation::Update(size_t index_i, Real dt)
		{
			rho_n_[index_i] += drho_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		PressureRelaxationInnerOldroyd_B ::
			PressureRelaxationInnerOldroyd_B(BaseBodyRelationInner &inner_relation)
			: PressureRelaxationDissipativeRiemannInner(inner_relation),
			  tau_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_->base_particles_)->tau_),
			  dtau_dt_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_->base_particles_)->dtau_dt_) {}
		//=================================================================================================//
		void PressureRelaxationInnerOldroyd_B::Initialization(size_t index_i, Real dt)
		{
			PressureRelaxationDissipativeRiemannInner::Initialization(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void PressureRelaxationInnerOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxationDissipativeRiemannInner::Interaction(index_i, dt);

			Real rho_i = rho_n_[index_i];
			Matd tau_i = tau_[index_i];

			Vecd acceleration(0);
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				//elastic force
				acceleration += (tau_i + tau_[index_j]) * nablaW_ij * Vol_[index_j] / rho_i;
			}

			dvel_dt_[index_i] += acceleration;
		}
		//=================================================================================================//
		DensityRelaxationInnerOldroyd_B::
			DensityRelaxationInnerOldroyd_B(BaseBodyRelationInner &inner_relation)
			: DensityRelaxationDissipativeRiemannInner(inner_relation),
			  tau_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_->base_particles_)->tau_),
			  dtau_dt_(DynamicCast<ViscoelasticFluidParticles>(this, sph_body_->base_particles_)->dtau_dt_)
		{
			Oldroyd_B_Fluid *oldroy_b_fluid = DynamicCast<Oldroyd_B_Fluid>(this, sph_body_->base_material_);
			mu_p_ = oldroy_b_fluid->ReferencePolymericViscosity();
			lambda_ = oldroy_b_fluid->getReferenceRelaxationTime();
		}
		//=================================================================================================//
		void DensityRelaxationInnerOldroyd_B::Interaction(size_t index_i, Real dt)
		{
			DensityRelaxationDissipativeRiemannInner::Interaction(index_i, dt);

			Vecd vel_i = vel_n_[index_i];
			Matd tau_i = tau_[index_i];

			Matd stress_rate(0);
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];

				Matd velocity_gradient = -SimTK::outer((vel_i - vel_n_[index_j]), nablaW_ij) * Vol_[index_j];
				stress_rate += ~velocity_gradient * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
							   (~velocity_gradient + velocity_gradient) * mu_p_ / lambda_;
			}

			dtau_dt_[index_i] = stress_rate;
		}
		//=================================================================================================//
		void DensityRelaxationInnerOldroyd_B::Update(size_t index_i, Real dt)
		{
			DensityRelaxationDissipativeRiemannInner::Update(index_i, dt);

			tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		DensityRelaxationInnerShearThinning::
            DensityRelaxationInnerShearThinning(BaseBodyRelationInner &inner_relation)
            : DensityRelaxationDissipativeRiemannInner(inner_relation),
              shear_rate_(DynamicCast<ShearThinningFluidParticles>(this, sph_body_->base_particles_)->shear_rate_),
              mu_shear_(DynamicCast<ShearThinningFluidParticles>(this, sph_body_->base_particles_)->mu_shear_)
        {
            ShearThinningFluid *shear_thinning_fluid = DynamicCast<ShearThinningFluid>(this, sph_body_->base_material_);
            mu_ = material_->ReferenceViscosity();
            mu_0_ = shear_thinning_fluid->ReferenceZeroShearViscosity();
            lambda_ = shear_thinning_fluid->getReferenceOscillationTime();
            //a_ = shear_thinning_fluid->getReferenceParameterA();
            //b_ = shear_thinning_fluid->getReferenceParameterB();
        }
        //=================================================================================================//
        void DensityRelaxationInnerShearThinning::Interaction(size_t index_i, Real dt)
        {
            DensityRelaxationDissipativeRiemannInner::Interaction(index_i, dt);

            Vecd vel_i = vel_n_[index_i];
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];

            Matd shear_rate_tensor(0);
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                //viscous force
                Vecd nablaW_ij = inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
                Matd velocity_gradient_i = -SimTK::outer((vel_i - vel_n_[index_j]), nablaW_ij) * Vol_[index_j];
                shear_rate_tensor += 0.5 * (velocity_gradient_i + ~velocity_gradient_i);
            }
            shear_rate_[index_i] = sqrt(2.0) * shear_rate_tensor.norm();
			mu_shear_[index_i] = mu_ + (mu_0_ - mu_) * (1 + log(1 + lambda_ * shear_rate_[index_i])) / (1 + lambda_ * shear_rate_[index_i]);
            //mu_shear_[index_i] = mu_ + (mu_0_ - mu_) / std::pow(1 + std::pow(lambda_ * shear_rate_[index_i], b_), a_);
        }

	}
	//=================================================================================================//
}
//=================================================================================================//