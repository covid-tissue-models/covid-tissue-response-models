import sys

import math

import numpy as np
from cc3d.core.PySteppables import *
import random as rd
from statistics import mean
import os
from array import array
# Import project libraries and classes
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
# sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
# sys.path.append(os.environ["ViralInfectionVTM"])
# Import project libraries and classes
sys.path.append(os.path.dirname(__file__))
from Simulation.ViralInfectionVTMSteppableBasePy import *
import ViralInfectionVTMLib
# from ViralInfectionVTMModelInputs import *
from BatchRun import BatchRunLib

# Import toolkit
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from nCoVToolkit import nCoVUtils

from UniCellModelInputs import *


a = alpha
b = beta
# density = .1


# rs = 1

class UniCell_Explicit_Force_2DSteppable(ViralInfectionVTMSteppableBasePy):

    def __init__(self, frequency=1):

        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        # import ViralInfectionVTMModelInputs as ViralInfectionVTMModelInputs
        # import Models.DrugDosingModel.DrugDosingInputs as DrugDosingInputs
        import UniCellModelInputs# as UniCellModelInputs
        BatchRunLib.apply_external_multipliers(__name__, UniCellModelInputs)
        self.alpha = alpha
        self.beta = beta
        # self.memory = 0.4 #must be less than 0.5

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """

        # to_seed = density * self.dim.x * self.dim.y * self.

        self.seed_cells(numero)

        self.plot_win = self.add_new_plot_window(title='Theta',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Theta', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)

        self.plot_win.add_plot("Theta", style='Lines', color='red', size=5)
        theta0 = np.pi / 4 + np.pi
        for cell in self.cell_list:
            theta = np.random.uniform(0, 2 * np.pi)
            Fx = np.cos(theta0)
            Fy = np.sin(theta0)
            Fz = 0
            cell.dict["Old_pos"] = [cell.xCOM, cell.yCOM, cell.zCOM]
            cell.dict["ExForce"] = [Fx, Fy, Fz]
            cell.dict["Scale"] = 100.0
            cell.dict["Theta"] = theta0

            # Make sure ExternalPotential plugin is loaded
            cell.lambdaVecX = cell.dict["Scale"] * cell.dict["ExForce"][0]
            cell.lambdaVecY = cell.dict["Scale"] * cell.dict["ExForce"][1]
            cell.lambdaVecZ = cell.dict["Scale"] * cell.dict["ExForce"][2]

    def seed_cells(self, N):

        for i in range(N):
            cell = self.new_cell(self.UNICELL)
            x = np.random.uniform(2, self.dim.x - 1)
            y = np.random.uniform(2, self.dim.y - 1)
            empty = bool(self.cell_field[x, y, 0])
            while not empty:
                x = np.random.uniform(2, self.dim.x - 1)
                y = np.random.uniform(2, self.dim.y - 1)
                empty = bool(self.cell_field[x, y, 0])
            self.cell_field[x, y, 0] = cell
            self.cell_field[x + 1, y, 0] = cell
            self.cell_field[x, y + 1, 0] = cell
            self.cell_field[x - 1, y, 0] = cell
            self.cell_field[x, y - 1, 0] = cell

            self.cell_field[x + 1, y - 1, 0] = cell
            self.cell_field[x + 1, y + 1, 0] = cell
            self.cell_field[x - 1, y + 1, 0] = cell
            self.cell_field[x - 1, y - 1, 0] = cell

    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation
        
        :param mcs: current Monte Carlo step
        """
        if mcs % 1 == 0:
            for cell in self.cell_list:

                theta = np.random.vonmises(0., 4.)
                Fx_noise = np.cos(theta)
                Fy_noise = np.sin(theta)

                Current_pos = [cell.xCOM, cell.yCOM, cell.zCOM]

                Vx = Current_pos[0] - cell.dict["Old_pos"][0]
                Vy = Current_pos[1] - cell.dict["Old_pos"][1]
                Vz = Current_pos[2] - cell.dict["Old_pos"][2]

                Norm = np.sqrt(Vx * Vx + Vy * Vy + Vz * Vz)

                if Norm > 0:

                    theta_v = np.arctan2(Vy, Vx)
                    theta_f = np.arctan2(cell.lambdaVecY, cell.lambdaVecX)

                    dif = theta_f - theta_v
                    if dif > np.pi: dif = -2 * np.pi + dif
                    if dif < -np.pi: dif = 2 * np.pi + dif

                    theta_f += (1 - self.alpha) * dif
                    if theta_f > np.pi: theta_f = -2 * np.pi + theta_f
                    if theta_f < -np.pi: theta_f = 2 * np.pi + theta_f

                    theta_f += theta * self.beta
                    if theta_f > np.pi: theta_f = -2 * np.pi + theta_f
                    if theta_f < -np.pi: theta_f = 2 * np.pi + theta_f

                    cell.dict["Theta"] = theta_f

                    Fx = np.cos(theta_f)
                    Fy = np.sin(theta_f)
                    Fz = 0.

                    # Fx = cell.dict["ExForce"][0]*self.alfa + (1-self.alfa)*Vx/Norm
                    # Fy = cell.dict["ExForce"][1]*self.alfa + (1-self.alfa)*Vy/Norm
                    # Fz = cell.dict["ExForce"][2]*self.alfa + (1-self.alfa)*Vz/Norm

                    # Fx = Fx*self.beta + (1-self.beta)*Fx_noise
                    # Fy = Fy*self.beta + (1-self.beta)*Fy_noise

                    # self.plot_win.add_data_point("Theta", mcs, np.arctan2(Fy,Fx))

                    FNorm = np.sqrt(Fx * Fx + Fy * Fy + Fz * Fz)

                    # cell.dict["ExForce"] = [Fx/FNorm, Fy/FNorm, Fz/FNorm]
                    cell.dict["ExForce"] = [Fx, Fy, Fz]

                    cell.lambdaVecX = cell.dict["Scale"] * cell.dict["ExForce"][0]
                    cell.lambdaVecY = cell.dict["Scale"] * cell.dict["ExForce"][1]
                    # cell.lambdaVecZ = cell.dict["Scale"]*cell.dict["ExForce"][2]

                    # print(cell.dict["Old_pos"][:])
                    cell.dict["Old_pos"][0] = Current_pos[0]
                    cell.dict["Old_pos"][1] = Current_pos[1]

    def finish(self):
        """
        Called after the last MCS to wrap up the simulation
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """


class CalculationsSteppable(ViralInfectionVTMSteppableBasePy):
    def __init__(self, frequency=1):
        '''
        constructor
        '''
        SteppableBasePy.__init__(self, frequency)
        # PLACE YOUR CODE BELOW THIS LINE

    def start(self):
        self.Ncell = 0
        for compartments in self.clusters:
            self.Ncell += 1
        # self.file1 = open(r"D:\CompuCell3D-py3-64bit\Simulations\UniCell_Explicit_Force_2D\Output.txt", "a")
        self.output_path = Path(self.output_dir).joinpath("Output.txt")
        self.file1 = open(self.output_path, "a")

        self.file1.write("Distance \t p.p \n")

    def step(self, mcs):

        # calculate polarity and center of mass of a random cell
        dimx = self.dim.x
        dimy = self.dim.y
        cell_id = rd.randint(1, self.Ncell)  # CHOOSE A RANDOM CELL
        _cell = self.fetch_cell_by_id(cell_id)
        _vol = _cell.volume
        _cm = [_cell.xCOM, _cell.yCOM]

        _Polarity = [_cell.dict["ExForce"][0], _cell.dict["ExForce"][1]]

        # print(_cm)

        # print the polarity product and the distance between this cell and all the others
        for cell in self.cell_list:
            cell_vol = cell.volume
            cell_cm = [cell.xCOM, cell.yCOM]

            Polarity = [cell.dict["ExForce"][0], cell.dict["ExForce"][1]]
            CM = [cell.xCOM, cell.yCOM]

            if CM[0] - _cm[0] > dimx / 2.:
                correctx = -1
            elif CM[0] - _cm[0] < -dimx / 2.:
                correctx = 1
            else:
                correctx = 0
            if CM[1] - _cm[1] > dimy / 2.:
                correcty = -1
            elif CM[1] - _cm[1] < -dimy / 2.:
                correcty = 1
            else:
                correcty = 0

            CM = np.add(CM, [dimx * correctx, dimy * correcty])

            P_product = np.inner(Polarity, _Polarity)

            distance = np.linalg.norm(np.subtract(CM, _cm))

            if (100 > distance > 0):
                self.file1.write(str(distance) + "\t" + str(P_product) + "\n")
        if not mcs % 1000:
            self.file1.flush()
            os.fsync(self.file1)

    def finish(self):
        '''
        this function may be called at the end of simulation - used very infrequently though
        '''
        # PLACE YOUR CODE BELOW THIS LINE
        self.file1.flush()
        os.fsync(self.file1)
        self.file1.close()
        return

    def on_stop(self):
        '''
        this gets called each time user stops simulation
        '''
        # PLACE YOUR CODE BELOW THIS LINE

        return


class PersistentNeighborsSteppable(ViralInfectionVTMSteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):

        out_dir_name = "UCEF2D"
        print(out_dir_name)
        if not os.path.exists(out_dir_name): os.makedirs(out_dir_name)
        file_name = "PN_a" + str(alpha) + "_b" + str(beta)  # +"_rs"+str(rs)
        # self.output_path = str(Path(out_dir_name + "\\" + file_name))
        self.output_path = Path(self.output_dir).joinpath(file_name)
        self.file4 = open(self.output_path, 'a')
        # self.file4.write("DeltaT \t PN \n")

        self.samples = 100
        self.DTmin = 100

        for cell in self.cell_list:
            cell.dict["ListN"] = np.zeros((self.samples, 20))

        self.count1 = 0

    def step(self, mcs):
        waiting_time = 10000
        if mcs > waiting_time:
            if (mcs - waiting_time) % self.DTmin == 0:
                for cell in self.cell_list:
                    self.count2 = 0
                    for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                        if neighbor:
                            cell.dict["ListN"][self.count1][self.count2] = neighbor.id + 1
                            self.count2 += 1
                        else:
                            cell.dict["ListN"][self.count1][self.count2] = 1
                            self.count2 += 1
                self.count1 += 1

            if (mcs - waiting_time) % (self.samples * self.DTmin) == 0:

                for cell in self.cell_list:
                    List0 = cell.dict["ListN"][0][:]
                    List0 = [i for i in List0 if i != 0]
                    for count3 in range(self.samples):
                        List_dt = cell.dict["ListN"][count3][:]
                        dt = self.DTmin*count3
                        # CN_list = [x for x in np.concatenate((List0, List_dt)) if x not in List0 or x not in List_dt]
                        CN_list = [x for x in List0 if x not in List_dt]
                        # CN = len(CN_list) * 0.5
                        CN = len(CN_list)/len(List0)
                        # print(List0, List_dt)
                        # print(CN_list)
                        if not 1 in np.concatenate((List0,List_dt)):
                            self.file4.write(str(dt)+"\t"+str(CN)+"\n")
                        # self.file4.write(str(dt) + " " + str(CN) + "\n")

                self.count1 = 0
                for cell in self.cell_list:
                    cell.dict["ListN"] = np.zeros((self.samples, 10))
            if not mcs % 1000:
                self.file4.flush()
                os.fsync(self.file4)
    def finish(self):
        self.file4.flush()
        os.fsync(self.file4)
        self.file4.close()
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        return


class CollectivityCalcSteppable(ViralInfectionVTMSteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):

        self.Collectivity = []
        self.phi = []
        self.gama = []

    def step(self, mcs):
        if mcs > 10000 and mcs % 10 == 0:
            count = 0
            phi_x = 0
            phi_y = 0
            L = 0
            for cell in self.cell_list:
                count += 1
                Fx = np.cos(cell.dict["Theta"])
                Fy = np.sin(cell.dict["Theta"])
                # print (Fx, cell.dict["ExForce"][0])
                X = cell.xCOM
                Y = cell.yCOM
                L += (- Fx * (Y - self.dim.y / 2.) + Fy * (X - self.dim.x / 2.)) / np.sqrt(
                    (X - self.dim.x / 2.) ** 2 + (Y - self.dim.y / 2.) ** 2)
                phi_x += Fx
                phi_y += Fy
            phi_x /= count
            phi_y /= count
            P = np.sqrt(phi_x ** 2 + phi_y ** 2)
            L = np.sqrt(L ** 2) / count
            self.phi.append(P)
            self.gama.append(L)
            self.Collectivity.append(np.sqrt(P ** 2 + L ** 2))

    def finish(self):

        out_dir_name = "UCEF2D"
        print(out_dir_name)
        if not os.path.exists(out_dir_name): os.makedirs(out_dir_name)
        file_name = "a_phi_gama_col_a" + str(a) + "_b" + str(b)  # +"_rs"+str(rs)
        # self.output_path = str(Path(out_dir_name + "\\" + file_name))
        self.output_path = Path(self.output_dir).joinpath(file_name)
        self.file3 = open(self.output_path, 'a')
        self.file3.write(str(a) + "\t" + str(mean(self.phi)) + "\t" + str(mean(self.gama)) + "\t" + str(
            mean(self.Collectivity)) + "\n")
        # self.file3.write("\t"+str(mean(self.gama))+"\n")
        # self.file3.write("\t"+str(mean(self.Collectivity))+"\n")
        self.file3.flush()
        os.fsync(self.file3)
        self.file3.close()

        return

    def on_stop(self):

        return


class Position_OutputSteppable(ViralInfectionVTMSteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):

        out_dir_name = "UCEF2D"
        print(out_dir_name)
        if not os.path.exists(out_dir_name): os.makedirs(out_dir_name)
        file_name = "_a" + str(a) + "_b" + str(b)  # +"_rs"+str(rs)
        # self.output_path = str(Path(out_dir_name + "\\" + file_name))
        self.output_path = Path(self.output_dir).joinpath(file_name)

        # self.file4 = open(self.output_path, 'wb')
        # self.file4.write("MCS \t".encode())
        for cell in self.cell_list:
            # self.file4.write(f"X {cell.id}\t Y {cell.id}\t Px {cell.id} \t Py {cell.id}\t".encode())
                # "X " + str(cell.id) + "\t Y " + str(cell.id) + "\t Px " + str(cell.id) + "\t Py " + str(cell.id) + "\t")
            cell.dict["cx"] = 0
            cell.dict["cy"] = 0
            cell.dict["Old_pos2"] = [cell.xCOM, cell.yCOM, cell.zCOM]
        # self.file4.write("\n".encode())
        fname = f"medium_contact_a{a}_b{b}.dat"
        self.contact_file = open(Path(self.output_dir).joinpath(fname), "w+")
        self.contact_file.write("mcs, area\n")

    def step(self, mcs):

        if mcs > 100:
            # self.file4.write(f"{mcs}\t".encode())
            medium_contact = 0
            llist = []
            for cell in self.cell_list:
                neighbor_list = self.get_cell_neighbor_data_list(cell)
                medium_contact += neighbor_list.common_surface_area_with_cell_types(cell_type_list=[0])
                
                current_pos = [cell.xCOM, cell.yCOM, cell.zCOM]
                if cell.xCOM - cell.dict["Old_pos2"][0] > self.dim.x * 0.5: cell.dict["cx"] -= 1
                if cell.xCOM - cell.dict["Old_pos2"][0] < -self.dim.x * 0.5: cell.dict["cx"] += 1
                if cell.yCOM - cell.dict["Old_pos2"][1] > self.dim.y * 0.5: cell.dict["cy"] -= 1
                if cell.yCOM - cell.dict["Old_pos2"][1] < -self.dim.y * 0.5: cell.dict["cy"] += 1
                # self.file4.write(f'{cell.xCOM + self.dim.x * cell.dict["cx"]}\t{cell.yCOM + self.dim.y * cell.dict["cy"]}\t'.encode())
                # str(cell.xCOM + self.dim.x * cell.dict["cx"]) + "\t" + str(
                    # cell.yCOM + self.dim.y * cell.dict["cy"]) + "\t")
                # self.file4.write(f"{cell.lambdaVecX}\t{cell.lambdaVecY}\t".encode())
                # str(cell.lambdaVecX) + "\t" + str(cell.lambdaVecY) + "\t")
                cell.dict["Old_pos2"][:] = current_pos[:]
                # llist.append(cell.xCOM+self.dim.x*cell.dict["cx"])
                # llist.append(cell.yCOM+self.dim.y*cell.dict["cy"])
                # llist.append(cell.lambdaVecX)
                # llist.append(cell.lambdaVecY)
                
            self.contact_file.write(f"{mcs}, {medium_contact}\n")
            # arr = array("d",llist)
            # arr.tofile(self.file4)
            #self.file4.write("\n")
            if not mcs % 1000:
                # self.file4.flush()
                # os.fsync(self.file4)
                self.contact_file.flush()
                os.fsync(self.contact_file)
            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        # self.file4.flush()
        # os.fsync(self.file4)
        # self.file4.close()
        self.contact_file.flush()
        os.fsync(self.contact_file)
        self.contact_file.close()
        return

    def on_stop(self):
        # this gets called each time user stops simulation
        
        return
