/**
 * @file   3d_elastic_valve.cpp
 * @brief  3D simplified aortic heart valve against a pulsatile flow 
 * @author Massoud Rezavand, Virtonomy GmbH
 */
#include "sphinxsys.h"
using namespace SPH;

Real scale = 0.001;
Real DH = 6.2*scale;
Real DL = 10 * DH/2.0; 
Real DL1 = DH;	
Real DL2 = 0.0;
int number_of_particles = 40;		    		     		
Real resolution_ref = DH/number_of_particles; 			/**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0;	        /**< Sponge region to impose inflow condition. */
//Real DL_sponge = DL1;
Real BW = resolution_ref * 4.0; 			    /**< Boundary width. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-0.5 * DH - BW, -DL1 - BW, -0.5 * DH -BW),
	Vec3d(0.5 * DH + BW, DL + DL2 + BW, 0.5 * DH + BW));
BoundingBox fluid_body_domain_bounds(Vec3d(-0.5 * DH, -DL1, -0.5 * DH),
	Vec3d(0.5 * DH, DL + DL2, 0.5 * DH));
Vec3d emitter_halfsize = Vec3d(0.5 * DH, BW/2., 0.5 * DH);
Vec3d emitter_translation = Vec3d(0.0, -(DL1 - BW/2.), 0.0);
Vec3d inlet_buffer_halfsize = Vec3d(0.5 * DH, 0.5 * DL_sponge, 0.5 * DH);
Vec3d inlet_buffer_translation = Vec3d(0.0, -(DL1 - DL_sponge/2.), 0.0);
/**
 * @brief Material properties of blood as the fluid part.
 */
Real rho0_f = 1050.0;	     	            /**< Density of Blood. */
Real U_f = 2.e-6/(Pi/4*DH*DH);		        	            /**< Characteristic velocity. */
Real c_f = 10.0 * U_f * 4.;      	            /**< Speed of sound. */
Real mu_s = 3.6e-3;	                        /**< Approximate Dynamics viscosity of Blood. */
Real mu_p = 4.e-4;
Real mu_f = mu_s;
Real lambda_1 = 0.06;
Real lambda = 8.2;
Real mu_0 = 1.6e-1;
Real a = 1.23;
Real b = 0.64;			    

/** define the fluid body */
class WaterBlock : public ComplexShape
{
    public:
        WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
        {
            add<TriangleMeshShapeSTL>("./input/stenosed_vessel_blood_2R.stl", Vecd(0., 0., 0.)*scale, scale);  
        }
};

/** define the static solid wall boundary */
class WallBoundary : public ComplexShape
{
public:
    WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>("./input/stenosed_vessel_wall_40_2R.stl", Vecd(0., 0., 0.)*scale, scale);
    }
};

/**
 * @brief define velocity buffer
 */
class ParabolicInflow : public fluid_dynamics::InflowBoundaryCondition
{
    Real u_ave_, u_ref_, t_ref_;
public:
	ParabolicInflow(FluidBody &body, BodyAlignedBoxByCell &aligned_box_part)
		: InflowBoundaryCondition(body, aligned_box_part)
        ,u_ave_(0), u_ref_(U_f), t_ref_(1.0)
	{}
	Vecd getTargetVelocity(Vecd& position, Vecd& velocity)
	{
		Real u = velocity[0];
		Real v = velocity[1];
        Real w = velocity[2];
		if (position[1] < halfsize_[1]) {
            Real a = sqrt(position[0]*position[0]+position[2]*position[2]);
			u = 0.0;
            v = 2. * U_f * (1 - a * a / halfsize_[0] / halfsize_[0]);
            //v = 2. * u_ave_ * (1 - a * a / halfsize_[0] / halfsize_[0]);
            w = 0.0;
		}
		return Vecd(u, v, w);
	}

    //void setupDynamics(Real dt = 0.0) override
    //{
    //    Real run_time = GlobalStaticVariables::physical_time_;
    //    u_ave_ = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
    //}
};

/** fluid observer particle generator */
class FluidObserverParticleGenerator1 : public ObserverParticleGenerator
{
public:
    FluidObserverParticleGenerator1(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        int ny = 21;
        int nz = 21;

        vector<double> radius(ny, DH/2);
        radius[5] = 2.325 * scale;
        radius[6] = DH / 4;
        radius[7] = 1.785 * scale;
        radius[8] = 2.325 * scale;
        radius[9] = 2.865 * scale;

        for (int i = 0; i < ny; i++)
        {
            Real y = DL / (ny - 1) * i;
            for (int j = 0; j < nz; j++)
            {
                Real z = -radius[i] + 2 * radius[i] / (nz - 1) * j;
                positions_.push_back(Vecd(0.0, y, z));
            }
        }
    }
};

class FluidObserverParticleGenerator2 : public ObserverParticleGenerator
{
public:
    FluidObserverParticleGenerator2(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        int ny = 101;
        int nz = 21;
        for (int i = 0; i < ny; i++)
        {
            Real y = DL / (ny - 1) * i;
            for (int j = 0; j < nz; j++)
            {
                Real z = -DH/2. + 2 * DH/2. / (nz - 1) * j;
                positions_.push_back(Vecd(0.0, y, z));
            }
        }
    }
};


/**
 * @brief 	Main program starts here.
 */
int main()
{
    /**
     * @brief Build up -- a SPHSystem --
     */
    SPHSystem system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation**/
    system.run_particle_relaxation_ = false;
    /** Tag for restarting with relaxed body-fitted particles**/
    system.reload_particles_ = true;
    /** Output.*/
    InOutput in_output(system);
    /** Set the starting time. */
    GlobalStaticVariables::physical_time_ = 0.0;
    /** Tag for computation from restart files. 0: not from restart files. */
    system.restart_step_ = 0;
    /**
     * @brief Material property, partilces and body creation of fluid.
     */
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.setBodyDomainBounds(fluid_body_domain_bounds);
    //water_block.defineParticlesAndMaterial<ViscoelasticFluidParticles, Oldroyd_B_Fluid>(rho0_f, c_f, mu_s, lambda_1, mu_p);
	water_block.defineParticlesAndMaterial<ShearThinningFluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    //water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();
    /**
     * @brief 	Particle and body creation of wall boundary.
     */
    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    //wall_boundary.generateParticles<ParticleGeneratorLattice>();

    if (!system.run_particle_relaxation_ && system.reload_particles_)
    {
        wall_boundary.generateParticles<ParticleGeneratorReload>(in_output, wall_boundary.getBodyName());
    }
    else
    {
        wall_boundary.defineBodyLevelSetShape()->writeLevelSet(wall_boundary);
        wall_boundary.generateParticles<ParticleGeneratorLattice>();
    }
    /**
     * @brief 	Particle and body creation of fluid observer.
     */
    ObserverBody fluid_observer1(system, "FluidObserver1");
    fluid_observer1.generateParticles<FluidObserverParticleGenerator1>();
    
    ObserverBody fluid_observer2(system, "FluidObserver2");
    fluid_observer2.generateParticles<FluidObserverParticleGenerator2>();

    /** Output the body states. */
    BodyStatesRecordingToVtp body_states_recording(in_output, system.real_bodies_);
    /** Output the body states for restart simulation. */
    RestartIO		restart_io(in_output, system.real_bodies_);
    /** topology */
	BodyRelationInner water_body_inner(water_block);
    ComplexBodyRelation water_block_complex(water_body_inner, {&wall_boundary});
    BodyRelationContact fluid_observer_contact1(fluid_observer1, {&water_block});
    BodyRelationContact fluid_observer_contact2(fluid_observer2, {&water_block});

    //----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	if (system.run_particle_relaxation_)
    {
        BodyRelationInner wall_bound_model_inner(wall_boundary);
	    RandomizeParticlePosition random_wall_particles(wall_boundary);

	    BodyStatesRecordingToVtp write_relaxed_particles(in_output, system.real_bodies_);
	    ReloadParticleIO write_particle_reload_files(in_output, { &wall_boundary});

	    relax_dynamics::RelaxationStepInner relaxation_step_inner_wall(wall_bound_model_inner, true);
        //----------------------------------------------------------------------
	    //	Particle relaxation starts here.
	    //----------------------------------------------------------------------
	    random_wall_particles.parallel_exec(0.25);
	    relaxation_step_inner_wall.surface_bounding_.parallel_exec();
	    write_relaxed_particles.writeToFile(0);
        wall_boundary.updateCellLinkedList();
	    //----------------------------------------------------------------------
	    //	Particle relaxation time stepping start here.
	    //----------------------------------------------------------------------
	    int ite_p = 0;
	    while (ite_p < 1000)
	    {
            relaxation_step_inner_wall.parallel_exec();
	    	ite_p += 1;
	    	if (ite_p % 100 == 0)
	    	{
	    		std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the model N = " << ite_p << "\n";
                write_relaxed_particles.writeToFile(ite_p);
	    	}
	    }
	    write_particle_reload_files.writeToFile(0);
        return 0;
    }
        /**
     * @brief 	Methods used for time stepping.
     */
     /** Initialize particle acceleration. */
    TimeStepInitialization 	initialize_a_fluid_step(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
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
	/** Evaluation of density by freestream approach. */
	fluid_dynamics::DensitySummationFreeStreamComplex update_density_by_summation(water_block_complex);
	/** We can output a method-specific particle data for debug */
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<int>("SurfaceIndicator");

    //water_block.addBodyStateForRecording<Matd>("ElasticStress");
    //water_block.addBodyStateForRecording<Matd>("ElasticStressChangeRate");
    /**
     * @brief 	Algorithms of fluid dynamics.
     */
     /** Evaluation of density by summation approach. */
    //fluid_dynamics::DensitySummationComplex		update_density_by_summation(water_block_complex);
    /** Time step size without considering sound wave speed. */
    fluid_dynamics::AdvectionTimeStepSize 	get_fluid_advection_time_step_size(water_block, U_f);
    /** Time step size with considering sound wave speed. */
    fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
    //ParabolicInitialVelocity initial_condition(water_block);
    /** Pressure relaxation algorithm without Riemann solver for viscous flows. */
    fluid_dynamics::PressureRelaxationRiemannWithWall	pressure_relaxation(water_block_complex);
    //fluid_dynamics::PressureRelaxationWithWallOldroyd_B pressure_relaxation(water_block_complex);
    /** Pressure relaxation algorithm by using position verlet time stepping. */
    fluid_dynamics::DensityRelaxationRiemannWithWall	density_relaxation(water_block_complex);
    //fluid_dynamics::DensityRelaxationWithWallOldroyd_B density_relaxation(water_block_complex);
    /** Computing viscous acceleration. */
    fluid_dynamics::UpdateViscosity update_viscosity(water_body_inner, mu_0, lambda, a, b);
    water_block.addBodyStateForRecording<Real>("ShearRate");
    water_block.addBodyStateForRecording<Real>("Viscosity");
    fluid_dynamics::ShearThinningViscousAccelerationWithWall	viscous_acceleration(water_block_complex);
    //fluid_dynamics::ViscousAccelerationWithWall	viscous_acceleration(water_block_complex);
    /** Impose transport velocity. */
    fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
    /** Define outlet face to */
    OpenBoundaryConditionInAxisDirection transfer_to_buffer_particles(water_block, yAxis, positiveDirection);

    /** Define the methods for observations of the simulation. */
    ObservedQuantityRecording<Vecd> write_fluid_velocity1("Velocity", in_output, fluid_observer_contact1);
    ObservedQuantityRecording<Vecd> write_fluid_velocity2("Velocity", in_output, fluid_observer_contact2);
    /**
     * @brief Setup geomtry and initial conditions.
     */
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    wall_boundary_normal_direction.parallel_exec();
    //initial_condition.exec();
    /**
     * @brief The time stepping starts here.
     */
     /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
    if (system.restart_step_ != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
        water_block.updateCellLinkedList();
        water_block_complex.updateConfiguration();
    }
    /** Output the starting states of bodies. */
    body_states_recording.writeToFile(0);
    //write_states.writeToFile(0);
    /**
     * @brief 	Basic parameters.
     */
    size_t number_of_iterations = system.restart_step_;
    int screen_output_interval = 100;
    int restart_output_interval = screen_output_interval*10;
    Real End_Time 		= 0.7; 	/**< End time. */
    Real D_Time = End_Time / 50.0;	/**< time stamps for output. */
    Real Dt = 0.0;			        /**< Default advection time step sizes. */
    Real dt = 0.0; 			        /**< Default acoustic time step sizes. */
    /** statistics for computing CPU time. */
    tick_count t1 = tick_count::now();
    tick_count::interval_t interval;
    /**
     * @brief 	Main loop starts here.
     */
    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < D_Time)
        {
            /** Acceleration due to viscous force and correct with TV. */
            initialize_a_fluid_step.parallel_exec();
            Dt = get_fluid_advection_time_step_size.parallel_exec();
            inlet_outlet_surface_particle_indicator.parallel_exec();
            update_density_by_summation.parallel_exec();
            update_viscosity.parallel_exec();
            viscous_acceleration_and_transport_correction.parallel_exec(Dt);

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt - relaxation_time);
                pressure_relaxation.parallel_exec(dt);
                emitter_buffer_inflow_condition.exec();
                density_relaxation.parallel_exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0 && number_of_iterations != system.restart_step_)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

            emitter_inflow_injecting.exec();
            transfer_to_buffer_particles.particle_type_transfer.parallel_exec();
            water_block.updateCellLinkedList();
            water_block_complex.updateConfiguration();
        }

        tick_count t2 = tick_count::now();
        body_states_recording.writeToFile();
        fluid_observer_contact1.updateConfiguration();
        fluid_observer_contact2.updateConfiguration();
        write_fluid_velocity1.writeToFile(number_of_iterations);
        write_fluid_velocity2.writeToFile(number_of_iterations);
        tick_count t3 = tick_count::now();
        interval += t3 - t2;
    }
    tick_count t4 = tick_count::now();

    tick_count::interval_t tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall clock time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
