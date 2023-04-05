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

constexpr double cPi = 3.141592653589793;
constexpr double c2PiO3 = 2.0*cPi/3.0;

int get_index(int i_turbine, int n_turbines, int i_blade, int n_blades, int i_node, int n_nodes)
{
    return i_node + n_nodes*(i_blade + n_blades*i_turbine);
}

template<class dtype>
dtype* convert_to_pointer(py::array_t<dtype> array)
{
    return static_cast<dtype *>(array.request().ptr);
}

void set_velocities(fast::OpenFAST& self,
    py::array_t<double> nacelle_velocities_arr, 
    py::array_t<double> blade_velocities_x_arr, py::array_t<double> blade_velocities_y_arr, py::array_t<double> blade_velocities_z_arr,
    py::array_t<double> tower_velocities_x_arr, py::array_t<double> tower_velocities_y_arr, py::array_t<double> tower_velocities_z_arr)
{


    auto nacelle_velocities = convert_to_pointer(nacelle_velocities_arr);
    auto blade_velocities_x = convert_to_pointer(blade_velocities_x_arr);
    auto blade_velocities_y = convert_to_pointer(blade_velocities_y_arr);
    auto blade_velocities_z = convert_to_pointer(blade_velocities_z_arr);
    
    auto tower_velocities_x = convert_to_pointer(tower_velocities_x_arr);
    auto tower_velocities_y = convert_to_pointer(tower_velocities_y_arr);
    auto tower_velocities_z = convert_to_pointer(tower_velocities_z_arr);

    //for now assume same nBlades and nBladeNodes at all turbines
    int nTurbines = self.get_nTurbinesProc();
    int nBlades = self.get_numBlades(0);
    int nBladeNodes = self.get_numForcePtsBlade(0);
    int nTowerNodes = self.get_numForcePtsTwr(0);

    int iNode = 0;
    for(int i_turbine=0; i_turbine<nTurbines; i_turbine++)
    {
        int iTurbGlob = self.get_globalTurbNo(i_turbine);

        // Nacelle
        double vel[3] = {nacelle_velocities[i_turbine*3], nacelle_velocities[i_turbine*3+1], nacelle_velocities[i_turbine*3+2]};
        self.setVelocityForceNode(vel, iNode, iTurbGlob);
        iNode++;

        // Blades
        for(int i_blade=0; i_blade<nBlades; i_blade++)
        {
            for(int i_blade_node=0; i_blade_node<nBladeNodes; i_blade_node++)
            {
                int index = get_index(i_turbine, nTurbines, i_blade, nBlades, i_blade_node, nBladeNodes);
                double vel[3] = {blade_velocities_x[index], blade_velocities_y[index], blade_velocities_z[index]};
                self.setVelocityForceNode(vel, iNode, iTurbGlob);
                iNode++;
            }
        }

        // Tower
        for(int i_tower_node=0; i_tower_node<nTowerNodes; i_tower_node++)
        {
            int index = get_index(i_turbine, nTurbines, 0,1, i_tower_node, nTowerNodes);
            double vel[3] = {tower_velocities_x[index], tower_velocities_y[index], tower_velocities_z[index]};
            self.setVelocityForceNode(vel, iNode, iTurbGlob);
            iNode++;
        }
    }
}

void get_forces(fast::OpenFAST& self,
    py::array_t<double> nacelle_forces_arr,
    py::array_t<double> blade_forces_x_arr, py::array_t<double> blade_forces_y_arr, py::array_t<double> blade_forces_z_arr,
    py::array_t<double> tower_forces_x_arr, py::array_t<double> tower_forces_y_arr, py::array_t<double> tower_forces_z_arr)
{

    auto nacelle_forces = convert_to_pointer(nacelle_forces_arr);
    auto blade_forces_x = convert_to_pointer(blade_forces_x_arr);
    auto blade_forces_y = convert_to_pointer(blade_forces_y_arr);
    auto blade_forces_z = convert_to_pointer(blade_forces_z_arr);
    
    auto tower_forces_x = convert_to_pointer(tower_forces_x_arr);
    auto tower_forces_y = convert_to_pointer(tower_forces_y_arr);
    auto tower_forces_z = convert_to_pointer(tower_forces_z_arr);

    //for now assume same nBlades and nBladeNodes at all turbines
    int nTurbines = self.get_nTurbinesProc();
    int nBlades = self.get_numBlades(0);
    int nBladeNodes = self.get_numForcePtsBlade(0);
    int nTowerNodes = self.get_numForcePtsTwr(0);


    int iNode = 0;
    for(int i_turbine=0; i_turbine<nTurbines; i_turbine++)
    {
        int iTurbGlob = self.get_globalTurbNo(i_turbine);

        double force[3];
        self.getForce(force, iNode, iTurbGlob);
        nacelle_forces[3*i_turbine  ] = -force[0];
        nacelle_forces[3*i_turbine+1] = -force[1];
        nacelle_forces[3*i_turbine+2] = -force[2];

        iNode++;

        for(int i_blade=0; i_blade<nBlades; i_blade++)
        {

            for(int i_blade_node=0; i_blade_node<nBladeNodes; i_blade_node++)
            {
                int index = get_index(i_turbine, nTurbines, i_blade, nBlades, i_blade_node, nBladeNodes);
                double forces[3];
                self.getForce(forces, iNode, iTurbGlob);
                iNode++;
                blade_forces_x[index] = -forces[0];
                blade_forces_y[index] = -forces[1];
                blade_forces_z[index] = -forces[2];
            }
        }

        for(int i_tower_node=0; i_tower_node<nTowerNodes; i_tower_node++)
        {
            int index = get_index(i_turbine, nTurbines, 0, 1, i_tower_node, nTowerNodes);
            double forces[3];
            self.getForce(forces, iNode, iTurbGlob);
            iNode++;
            tower_forces_x[index] = -forces[0];
            tower_forces_y[index] = -forces[1];
            tower_forces_z[index] = -forces[2];
        }
    }
}

void get_coordinates(fast::OpenFAST& self,
    py::array_t<double> nacelle_coordinates_arr,
    py::array_t<double> blade_coordinates_x_arr, py::array_t<double> blade_coordinates_y_arr, py::array_t<double> blade_coordinates_z_arr,
    py::array_t<double> tower_coordinates_x_arr, py::array_t<double> tower_coordinates_y_arr, py::array_t<double> tower_coordinates_z_arr)
{
    auto nacelle_coordinates = convert_to_pointer(nacelle_coordinates_arr);
    auto blade_coordinates_x = convert_to_pointer(blade_coordinates_x_arr);
    auto blade_coordinates_y = convert_to_pointer(blade_coordinates_y_arr);
    auto blade_coordinates_z = convert_to_pointer(blade_coordinates_z_arr);
    
    auto tower_coordinates_x = convert_to_pointer(tower_coordinates_x_arr);
    auto tower_coordinates_y = convert_to_pointer(tower_coordinates_y_arr);
    auto tower_coordinates_z = convert_to_pointer(tower_coordinates_z_arr);

    //for now assume same nBlades and nBladeNodes at all turbines
    int nTurbines = self.get_nTurbinesProc();
    int nBlades = self.get_numBlades(0);
    int nBladeNodes = self.get_numForcePtsBlade(0);
    int nTowerNodes = self.get_numForcePtsTwr(0);

    int iNode = 0;
    for(int i_turbine=0; i_turbine<nTurbines; i_turbine++)
    {
        int iTurbGlob = self.get_globalTurbNo(i_turbine);
        double hubPos[3];
        self.getHubPos(hubPos, iTurbGlob);
        double coordinates[3];
        self.getForceNodeCoordinates(coordinates, iNode, iTurbGlob);
        nacelle_coordinates[3*i_turbine  ] = coordinates[0] - hubPos[0];
        nacelle_coordinates[3*i_turbine+1] = coordinates[1] - hubPos[1];
        nacelle_coordinates[3*i_turbine+2] = coordinates[2] - hubPos[2];

        iNode++;

        for(int i_blade=0; i_blade<nBlades; i_blade++)
        {

            for(int i_blade_node=0; i_blade_node<nBladeNodes; i_blade_node++)
            {
                int index = get_index(i_turbine, nTurbines, i_blade, nBlades, i_blade_node, nBladeNodes);
                double coordinates[3];
                self.getForceNodeCoordinates(coordinates, iNode, iTurbGlob);
                iNode++;
                blade_coordinates_x[index] = coordinates[0] - hubPos[0];
                blade_coordinates_y[index] = coordinates[1] - hubPos[1];
                blade_coordinates_z[index] = coordinates[2] - hubPos[2];
            }
        }

        for(int i_tower_node=0; i_tower_node<nTowerNodes; i_tower_node++)
        {
            int index = get_index(i_turbine, nTurbines, 0, 1, i_tower_node, nTowerNodes);
            double coordinates[3];
            self.getForceNodeCoordinates(coordinates, iNode, iTurbGlob);
            iNode++;
            tower_coordinates_x[index] = coordinates[0];
            tower_coordinates_y[index] = coordinates[1];
            tower_coordinates_z[index] = coordinates[2];
        }
    }
}

void get_azimuths(fast::OpenFAST& self, py::array_t<double> azimuths_arr)
{

    double* azimuths = convert_to_pointer(azimuths_arr);
    int nTurbines = self.get_nTurbinesProc();
    int nBlades = self.get_numBlades(0);
    int nBladeNodes = self.get_numForcePtsBlade(0);
    int nTowerNodes = self.get_numForcePtsTwr(0);

    for(int i_turbine=0; i_turbine<nTurbines; i_turbine++)
    {
        int iTurbGlob = self.get_globalTurbNo(i_turbine);
        int index = get_index(i_turbine, nTurbines, 0, nBlades, 0, nBladeNodes);
        int iNode = i_turbine*(1+nBlades*nBladeNodes+nTowerNodes) + 2;

        double coords_blade_0[3];
        self.getForceNodeCoordinates(coords_blade_0, iNode, iTurbGlob);
        iNode = iNode + nBladeNodes;

        double coords_blade_1[3];
        self.getForceNodeCoordinates(coords_blade_1, iNode, iTurbGlob);
        iNode = iNode + nBladeNodes;

        double coords_blade_2[3];
        self.getForceNodeCoordinates(coords_blade_2, iNode, iTurbGlob);

        double center_x = (coords_blade_0[0] + coords_blade_1[0] + coords_blade_2[0]) / 3.0;
        double center_y = (coords_blade_0[1] + coords_blade_1[1] + coords_blade_2[1]) / 3.0;
        double center_z = (coords_blade_0[2] + coords_blade_1[2] + coords_blade_2[2]) / 3.0;

        double az0 = atan2( -(coords_blade_0[1]-center_y), coords_blade_0[2]-center_z);

        double az1 = atan2( -(coords_blade_1[1]-center_y), coords_blade_1[2]-center_z) - c2PiO3;
        az1 = az1 + (az1 < -cPi ? 2*cPi : 0);

        double az2 = atan2( -(coords_blade_2[1]-center_y), coords_blade_2[2]-center_z) - 2.0 * c2PiO3;
        az2 = az2 + (az2 < -cPi ? 2*cPi : 0);

        azimuths[i_turbine] = (az0+az1+az2)/3.0;
    }
}

void get_omegas(fast::OpenFAST& self, py::array_t<double> omegas_arr)
{

    double* omegas = convert_to_pointer(omegas_arr);
    int nTurbines = self.get_nTurbinesProc();

    for(int i_turbine=0; i_turbine<nTurbines; i_turbine++)
    {
        int iTurbGlob = self.get_globalTurbNo(i_turbine);
        omegas[i_turbine] = self.computeRotorSpeed(iTurbGlob);
    }
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
        .def("set_velocities", &set_velocities, py::arg("nacelle_velocities"), py::arg("blade_velocities_x"), py::arg("blade_velocities_y"), py::arg("blade_velocities_z"), py::arg("tower_velocities_x"), py::arg("tower_velocities_y"), py::arg("tower_velocities_z"))
        .def("get_forces", &get_forces, py::arg("nacelle_forces"), py::arg("blade_forces_x"), py::arg("blade_forces_y"), py::arg("blade_forces_z"), py::arg("tower_forces_x"), py::arg("tower_forces_y"), py::arg("tower_forces_z"))
        .def("get_coordinates", &get_coordinates, py::arg("nacelle_coordinates"), py::arg("blade_coordinates_x"), py::arg("blade_coordinates_y"), py::arg("blade_coordinates_z"), py::arg("tower_coordinates_x"), py::arg("tower_coordinates_y"), py::arg("tower_coordinates_z"))
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
        .def("get_azimuths", &get_azimuths, py::arg("azimuths"))
        .def("get_omegas", &get_omegas, py::arg("omegas"))
        .def("is_time_zero", &fast::OpenFAST::isTimeZero)
        .def("is_dry_run", &fast::OpenFAST::isDryRun);
}