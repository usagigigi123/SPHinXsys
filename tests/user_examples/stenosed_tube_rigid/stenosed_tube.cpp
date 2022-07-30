/**
* @file 	fsi2.h
* @brief 	This is the case file for the test of fluid - structure interaction.
* @details  We consider a flow - induced vibration of an elastic beam behind a cylinder in 2D.
* @author 	Xiangyu Hu, Chi Zhangand Luhui Han
*/

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.001;
Real DH = 5.*scale;							/**< Tube diameter. */
Real L = 1.5 * scale;
Real h = 0.646 * scale;
Real DL = 3.* DH + (sqrt(3) + 1.) * h + L + 10. * DH;		/**< Tube length. */
Real DH2 = 1.5 * DH;                                    /**< Largest diameter after deformation. */
int number_of_particles = 10;	
Real resolution_ref = DH/number_of_particles;				/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = resolution_ref * 4.0;			/**< Boundary width, determined by specific layer of boundary particles. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-0.5 * DH2 - BW, -DL_sponge - BW, -0.5 * DH2 - BW),
	Vec3d(0.5 * DH2 + BW, DL + BW, 0.5 * DH2 + BW));
BoundingBox fluid_body_domain_bounds(Vec3d(-0.5 * DH2, -DL_sponge, -0.5 * DH2),
	Vec3d(0.5 * DH2, DL, 0.5 * DH2));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 755.0;											 /**< Density. */
Real U_f = 2. * 6.9e-6 /(Pi/4*DH*DH) / (1 - 2*h/DH) / (1 - 2*h/DH);								    /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;		   /**< Speed of sound. */
Real mu_f = 1.43e-3;;                                          /**< Infinite viscosity. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
Vec3d emitter_halfsize = Vec3d(0.5 * DH, BW/2., 0.5 * DH);
Vec3d emitter_translation = Vec3d(0.0, -DL_sponge + BW/2., 0.0);
Vec3d inlet_buffer_halfsize = Vec3d(0.5 * DH, 0.5 * DL_sponge, 0.5 * DH);
Vec3d inlet_buffer_translation = Vec3d(0.0, - DL_sponge/2., 0.0);

/** resolution which controls the quality of polygonalmesh created by geometry system */
int resolution(20);
//----------------------------------------------------------------------
//	Define case dependent geometrices
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
    public:
        WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
        {
            add<TriangleMeshShapeSTL>("./input/stenosed_tube_block_3D_10D_10.stl", Vecd(0., 0., 0.)*scale, scale); 
        }
};
class WallBoundary : public ComplexShape
{
public:
    WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>("./input/stenosed_tube_rigid_wall_3D_10D_10.stl", Vecd(0., 0., 0.)*scale, scale);
    }
};
/** create base as constrain shape. */
class WallBase : public BodyPartByParticle
{
public:
	WallBase(SPHBody &body, const std::string &body_part_name) : BodyPartByParticle(body, body_part_name)
	{
		TaggingParticleMethod tagging_particle_method = std::bind(&WallBase::tagManually, this, _1);
		tagParticles(tagging_particle_method);
	}
	virtual ~WallBase(){};
private:
	void tagManually(size_t index_i)
	{
		if (base_particles_->pos_n_[index_i][1] <= 1.3 * BW || base_particles_->pos_n_[index_i][1] >= DL - BW)
		{
			body_part_particles_.push_back(index_i);
		}
	}
};

/** Case dependent inflow boundary condition. */
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
	Real u_ave_, Q_mid, Q_amp, t_ref;

public:
	ParabolicInflow(FluidBody &fluid_body, BodyAlignedBoxByCell &aligned_box_part)
		: InflowBoundaryCondition(fluid_body, aligned_box_part),
		  u_ave_(0), t_ref(0.345), Q_mid(4.3), Q_amp(2.6) {}
	Vecd getTargetVelocity(Vecd &position, Vecd &velocity) override
	{
		Real u = velocity[0];
		Real v = velocity[1];
		Real w = velocity[2];
		if (position[1] < halfsize_[1])
		{
			Real a = sqrt(position[0]*position[0]+position[2]*position[2]);
			u = 0.0;
			v = 2. * u_ave_ * (1 - a * a / halfsize_[0] / halfsize_[0]);
			w = 0.0;
		}
		return Vecd(u, v, w);
	}
	void setupDynamics(Real dt = 0.0) override
	{
		Real run_time = GlobalStaticVariables::physical_time_;
		//Real Q = Q_mid;
		//Real Q = run_time < t_ref ? Q_mid : Q_amp * sin(2 * Pi * run_time / t_ref) + Q_mid;
		Real Q = Q_amp * sin(2 * Pi * run_time / t_ref) + Q_mid;
		u_ave_ = Q * 1.e-6 / (Pi/4*DH*DH);
	}
};
/** fluid observer particle generator */
class FluidObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit FluidObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** A line of measuring points at the entrance of the channel. */
		size_t number_observation_points = 21;
		Real y_mid = 3. * DH + ((1+sqrt(3))*h + L) * 0.5;
		Real y1 = y_mid + DH;
		Real y2 = y_mid + 2.5 * DH;
		Real y3 = y_mid + 4.3 * DH;
		/** the measuring locations */
		for (size_t i = 0; i < number_observation_points; ++i)
		{
			Real z = DH * (Real)i / (Real)(number_observation_points - 1) - DH/2.;
			positions_.push_back(Vec3d(0.0, y1, z));
			positions_.push_back(Vec3d(0.0, y2, z));
			positions_.push_back(Vec3d(0.0, y3, z));
		}
	}
};

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	sph_system.run_particle_relaxation_ = false;
	/** Tag for computation start with relaxed body fitted particles distribution. */
	sph_system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: start with initial condition. */
	sph_system.restart_step_ = 0;
// handle command line arguments
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av);
#endif
	/** output environment. */
	InOutput in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_block.setBodyDomainBounds(fluid_body_domain_bounds);
	water_block.defineBodyLevelSetShape()->writeLevelSet(water_block);
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
	wall_boundary.defineBodyLevelSetShape()->writeLevelSet(wall_boundary);
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();

	if (!sph_system.run_particle_relaxation_ && sph_system.reload_particles_)
    {
        wall_boundary.generateParticles<ParticleGeneratorReload>(in_output, wall_boundary.getBodyName());
		water_block.generateParticles<ParticleGeneratorReload>(in_output, water_block.getBodyName());
    }
    else
    {
        wall_boundary.generateParticles<ParticleGeneratorLattice>();
		water_block.generateParticles<ParticleGeneratorLattice>();
    }

	//ObserverBody beam_observer(sph_system, "BeamObserver");
	//beam_observer.generateParticles<ObserverParticleGenerator>(beam_observation_location);
	ObserverBody fluid_observer(sph_system, "FluidObserver");
	fluid_observer.generateParticles<FluidObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner wall_boundary_inner(wall_boundary);
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
	BodyRelationContact wall_boundary_contact(wall_boundary, {&water_block});
	//BodyRelationContact beam_observer_contact(beam_observer, {&insert_body});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.run_particle_relaxation_)
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
	    RandomizeParticlePosition random_wall_particles(wall_boundary);
		RandomizeParticlePosition random_water_particles(water_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_relaxed_particles(in_output, sph_system.real_bodies_);
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, sph_system.real_bodies_);
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner_wall(wall_boundary_inner);
		relax_dynamics::RelaxationStepInner relaxation_step_inner_water(water_block_complex.inner_relation_);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_wall_particles.parallel_exec(0.25);
	    relaxation_step_inner_wall.surface_bounding_.parallel_exec();
		random_water_particles.parallel_exec(0.25);
	    relaxation_step_inner_water.surface_bounding_.parallel_exec();
	    write_relaxed_particles.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner_wall.parallel_exec();
			relaxation_step_inner_water.parallel_exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
				write_relaxed_particles.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_fluid_step(water_block);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_density_by_summation(water_block_complex);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** Pressure relaxation using verlet time stepping. */
	/** Here, we do not use Riemann solver for pressure as the flow is viscous. */
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInner compute_vorticity(water_block_complex.inner_relation_);
	/** Emitter. */
    BodyAlignedBoxByParticle emitter(
		water_block, makeShared<AlignedBoxShape>(Transform3d(Vec3d(emitter_translation)), emitter_halfsize));
	fluid_dynamics::EmitterInflowInjecting emitter_inflow_injecting(water_block, emitter, 10, 1, true);
	/** Emitter condition. */
	BodyAlignedBoxByCell emitter_buffer(
		water_block, makeShared<AlignedBoxShape>(Transform3d(Vec3d(inlet_buffer_translation)), inlet_buffer_halfsize));
	ParabolicInflow emitter_buffer_inflow_condition(water_block, emitter_buffer);
	/** time-space method to detect surface particles. */
	fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationComplex
		inlet_outlet_surface_particle_indicator(water_block_complex);
	/** Periodic BCs in x direction. */
	OpenBoundaryConditionInAxisDirection transfer_to_buffer_particles(water_block, yAxis, positiveDirection);
	/** We can output a method-specific particle data for debug */
	water_block.addBodyStateForRecording<Real>("Density");
	water_block.addBodyStateForRecording<Real>("Mass");
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<AngularVecd>("VorticityInner");
	water_block.addBodyStateForRecording<int>("SurfaceIndicator");

	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(in_output, sph_system.real_bodies_);
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	ObservedQuantityRecording<Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	/** initialize cell linked lists for all bodies. */
	sph_system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	sph_system.initializeSystemConfigurations();
	/** computing surface normal direction for the wall. */
	wall_boundary_normal_direction.parallel_exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block.updateCellLinkedList();
		/** one need update configuration after periodic condition. */
		water_block_complex.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 2.0 * 0.345;			/**< End time. */
	Real D_Time = End_Time / 200.0; /**< time stamps for output. */
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_real_body_states.writeToFile();
	//write_beam_tip_displacement.writeToFile(number_of_iterations);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			inlet_outlet_surface_particle_indicator.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec(Dt);

			size_t inner_ite_dt = 0;
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);
				/** Fluid pressure relaxation */
				pressure_relaxation.parallel_exec(dt);
				/** Buffer velocity */
				emitter_buffer_inflow_condition.exec();
				/** Fluid density relaxation */
				density_relaxation.parallel_exec(dt);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				inner_ite_dt++;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt 
						  << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != sph_system.restart_step_)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			emitter_inflow_injecting.exec();
            transfer_to_buffer_particles.particle_type_transfer.parallel_exec();
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
		}

		tick_count t2 = tick_count::now();
		/** write run-time observation into file */
		compute_vorticity.parallel_exec();
		write_real_body_states.writeToFile();
		fluid_observer_contact.updateConfiguration();
		write_fluid_velocity.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}

