#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "OpenFAST.H"
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"

// TODO get pitch
namespace py = pybind11;


template<class dtype>
dtype* convert_to_pointer(py::array_t<dtype> array)
{
    return static_cast<dtype *>(array.request().ptr);
}


void set_velocities(
    fast::OpenFAST& self, 
    py::array_t<double> nacelle_velocities, 
    py::array_t<double> blade_velocities_x, py::array_t<double> blade_velocities_y, py::array_t<double> blade_velocities_z, 
    py::array_t<double> tower_velocities_x, py::array_t<double> tower_velocities_y, py::array_t<double> tower_velocities_z,
    py::array_t<double> radii, py::array_t<double> heights) 
{
    int numberOfExternBladeNodes = radii.shape()[0];
    int numberOfExternTowerNodes = heights.shape()[0];
    self.setAllLocalVelocitiesFromBladeFrame( 
        convert_to_pointer(nacelle_velocities), 
        convert_to_pointer(blade_velocities_x), convert_to_pointer(blade_velocities_y), convert_to_pointer(blade_velocities_z), 
        convert_to_pointer(tower_velocities_x), convert_to_pointer(tower_velocities_y), convert_to_pointer(tower_velocities_z),
        convert_to_pointer(radii), numberOfExternBladeNodes, convert_to_pointer(heights), numberOfExternTowerNodes);
}


void get_forces(
    fast::OpenFAST& self, 
    py::array_t<double> nacelle_forces, 
    py::array_t<double> blade_forces_x, py::array_t<double> blade_forces_y, py::array_t<double> blade_forces_z, 
    py::array_t<double> tower_forces_x, py::array_t<double> tower_forces_y, py::array_t<double> tower_forces_z,
    py::array_t<double> radii) 
{
    int numberOfInputNodes = radii.shape()[0];

    self.getAllLocalForcesInBladeFrame( 
        convert_to_pointer(nacelle_forces), 
        convert_to_pointer(blade_forces_x), convert_to_pointer(blade_forces_y), convert_to_pointer(blade_forces_z), 
        convert_to_pointer(tower_forces_x), convert_to_pointer(tower_forces_y), convert_to_pointer(tower_forces_z),
        convert_to_pointer(radii), numberOfInputNodes);
}

void get_coordinates(
    fast::OpenFAST& self, 
    py::array_t<double> nacelle_coordinates, 
    py::array_t<double> blade_coordinates_x, py::array_t<double> blade_coordinates_y, py::array_t<double> blade_coordinates_z, 
    py::array_t<double> tower_coordinates_x, py::array_t<double> tower_coordinates_y, py::array_t<double> tower_coordinates_z,
    py::array_t<double> radii) 
{
    int numberOfInputNodes = radii.shape()[0];

    self.getAllLocalCoordinatesInBladeFrame( 
        convert_to_pointer(nacelle_coordinates), 
        convert_to_pointer(blade_coordinates_x), convert_to_pointer(blade_coordinates_y), convert_to_pointer(blade_coordinates_z), 
        convert_to_pointer(tower_coordinates_x), convert_to_pointer(tower_coordinates_y), convert_to_pointer(tower_coordinates_z),
        convert_to_pointer(radii), numberOfInputNodes);
}


PYBIND11_MODULE(bindings, m) {
    py::class_<fast::globTurbineDataType>(m, "GlobTurbineDataType")
    .def(py::init<  int,
                    std::string,
                    std::string,
                    std::vector<double>,
                    std::vector<double>,
                    int,
                    int,
                    float,
                    float,
                    float>(), py::arg("TurbID"), 
                                py::arg("FASTInputFileName"), 
                                py::arg("FASTRestartFileName"), 
                                py::arg("TurbineBasePos"), 
                                py::arg("TurbineHubPos"), 
                                py::arg("numForcePtsBlade"), 
                                py::arg("numForcePtsTwr"), 
                                py::arg("nacelle_cd")=0.0, 
                                py::arg("nacelle_area")=0.0, 
                                py::arg("air_density")=0.0);
    
    py::enum_<fast::ActuatorNodeType>(m, "ActuatorNodeType")
        .value("HUB", fast::ActuatorNodeType::HUB)
        .value("BLADE", fast::ActuatorNodeType::BLADE)
        .value("TOWER", fast::ActuatorNodeType::TOWER);

    py::enum_<fast::simStartType>(m, "SimStartType")
        .value("INIT", fast::simStartType::init)
        .value("TRUE_RESTART", fast::simStartType::trueRestart)
        .value("RESTART_DRIVER_INIT_FAST", fast::simStartType::restartDriverInitFAST);

    py::class_<fast::fastInputs>(m, "FastInputs")
        .def(py::init<>())
        .def(py::init([](
            int nTurbinesGlob,
            bool dryRun,
            bool debug,
            double tStart,
            fast::simStartType simStart,
            int nEveryCheckPoint,
            double tMax,
            double dtFAST,
            bool scStatus,
            std::string scLibFile,
            std::vector<fast::globTurbineDataType>  globTurbineData) {
                fast::fastInputs fi;
                fi.comm = MPI_COMM_WORLD; // Hack to get around pybind complaining about incomplete type
                fi.nTurbinesGlob = nTurbinesGlob;
                fi.dryRun = dryRun;
                fi.debug = debug;
                fi.tStart = tStart;
                fi.simStart = simStart;
                fi.nEveryCheckPoint = nEveryCheckPoint;
                fi.tMax = tMax;
                fi.dtFAST = dtFAST;
                fi.scStatus = scStatus;
                fi.scLibFile = scLibFile;
                fi.globTurbineData = globTurbineData;
                return fi;
            }),
                py::arg("nTurbinesGlob"),
                py::arg("dryRun"),
                py::arg("debug"),
                py::arg("tStart"),
                py::arg("simStart"),
                py::arg("nEveryCheckPoint"),
                py::arg("tMax"),
                py::arg("dtFAST"),
                py::arg("scStatus"),
                py::arg("scLibFile"),
                py::arg("globTurbineData")
        );


    py::class_<fast::OpenFAST>(m, "OpenFAST")
        .def(py::init<>())
        .def("set_inputs", &fast::OpenFAST::setInputs, py::arg("fast_inputs"))
        .def("init", &fast::OpenFAST::init)
        .def("solution0", &fast::OpenFAST::solution0)
        .def("step", &fast::OpenFAST::step)
        .def("step_no_write", &fast::OpenFAST::stepNoWrite)
        .def("end", &fast::OpenFAST::end)
        .def("calc_nacelle_force", &fast::OpenFAST::calc_nacelle_force, py::arg("u"), py::arg("v"), py::arg("w"), py::arg("cd"), py::arg("area"), py::arg("rho"), py::arg("fx"), py::arg("fy"), py::arg("fz"))
        .def("set_velocities", &set_velocities, py::arg("nacelle_velocities"), py::arg("blade_velocities_x"), py::arg("blade_velocities_y"), py::arg("blade_velocities_z"), py::arg("tower_velocities_x"), py::arg("tower_velocities_y"), py::arg("tower_velocities_z"), py::arg("radii"), py::arg("heights"))
        .def("get_forces", &get_forces, py::arg("nacelle_forces"), py::arg("blade_forces_x"), py::arg("blade_forces_y"), py::arg("blade_forces_z"), py::arg("tower_forces_x"), py::arg("tower_forces_y"), py::arg("tower_forces_z"), py::arg("radii"))
        .def("get_coordinates", &get_coordinates, py::arg("nacelle_coordinates"), py::arg("blade_coordinates_x"), py::arg("blade_coordinates_y"), py::arg("blade_coordinates_z"), py::arg("tower_coordinates_x"), py::arg("tower_coordinates_y"), py::arg("tower_coordinates_z"), py::arg("radii"))
        .def("interpolate_velocity_at_force_nodes_to_velocity_nodes", &fast::OpenFAST::interpolateVel_ForceToVelNodes)
        .def("allocate_turbines_to_procs_simple", &fast::OpenFAST::allocateTurbinesToProcsSimple)
        .def("compute_torque_thrust", [](fast::OpenFAST& self, py::array_t<double> torque, py::array_t<double> thrust)
            {
                auto torque_array = convert_to_pointer(torque);
                auto thrust_array = convert_to_pointer(thrust);
                std::vector<double> torque_buffer(3);
                std::vector<double> thrust_buffer(3);
                for(int i_turbine=0; i_turbine<self.get_nTurbinesProc(); i_turbine++)
                {
                    self.computeTorqueThrust(i_turbine, torque_buffer, thrust_buffer);
                    torque_array[3*i_turbine  ] = torque_buffer[0];
                    torque_array[3*i_turbine+1] = torque_buffer[1];
                    torque_array[3*i_turbine+2] = torque_buffer[2];

                    thrust_array[3*i_turbine  ] = thrust_buffer[0];
                    thrust_array[3*i_turbine+1] = thrust_buffer[1];
                    thrust_array[3*i_turbine+2] = thrust_buffer[2];
                }
            }, py::arg("torque"), py::arg("thrust"))
        .def("get_azimuths", [](fast::OpenFAST& self, py::array_t<double> azimuths){self.getAllLocalAzimuths(convert_to_pointer(azimuths));}, py::arg("azimuths"))
        .def("get_rotor_speeds", [](fast::OpenFAST& self, py::array_t<double> rotor_speeds){self.getAllLocalRotorSpeeds(convert_to_pointer(rotor_speeds));}, py::arg("rotor_speeds"))
        .def("is_time_zero", &fast::OpenFAST::isTimeZero)
        .def("is_dry_run", &fast::OpenFAST::isDryRun);
}