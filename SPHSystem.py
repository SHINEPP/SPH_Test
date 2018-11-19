# -*- coding: utf-8 -*-

import math

PI = 3.1415926
INF = 1e-12
BOUNDARY = 0.0001


class Particle:
    def __init__(self):
        self.identifier = -1
        self.position = [0.0, 0.0, 0.0]  # 位置
        self.velocity = [0.0, 0.0, 0.0]  # 速率
        self.acceleration = [0.0, 0.0, 0.0]  # 加速度
        self.ev = [0.0, 0.0, 0.0]
        self.density = 0.0  # 密度
        self.pressure = 0.0  # 压力
        self.surf_norm = 0.0
        self.next_index = -1


class SPHSystem:
    def __init__(self):
        print 'init'
        self.__max_particle = 30000
        self.__num_particle = 0
        self.__kernel = 0.08
        self.__mass = 0.02

        self.__world_size = [1.0, 1.0, 1.0]
        self.__cell_size = self.__kernel

        self.__grid_size = [0, 0, 0]
        self.__grid_size[0] = int(math.ceil(self.__world_size[0] / self.__cell_size))
        self.__grid_size[1] = int(math.ceil(self.__world_size[1] / self.__cell_size))
        self.__grid_size[2] = int(math.ceil(self.__world_size[2] / self.__cell_size))
        self.__total_cell = self.__grid_size[0] * self.__grid_size[1] * self.__grid_size[2]

        self.__gravity = [0.0, -6.8, 0.0]
        self.__wall_damping = -0.5
        self.__rest_density = 1000.0  # 初始密度
        self.__gas_constant = 1.0  # 气体系数
        self.__viscosity = 6.5
        self.__time_step = 0.003
        self.__surf_norm = 6.0
        self.__surf_coe = 0.1

        self.__poly6_value = 315.0 / (64.0 * PI * math.pow(self.__kernel, 9))
        self.__spiky_value = -45.0 / (PI * math.pow(self.__kernel, 6))
        self.__viscosity_value = 45.0 / (PI * pow(self.__kernel, 6))

        self.__grad_poly6 = -945.0 / (32.0 * PI * math.pow(self.__kernel, 9))
        self.__lplc_poly6 = -945.0 / (8.0 * PI * math.pow(self.__kernel, 9))

        kernel_2 = self.__kernel * self.__kernel
        self.__self_density = self.__mass * self.__poly6_value * math.pow(self.__kernel, 6)
        self.__self_lplc_color = self.__lplc_poly6 * self.__mass * kernel_2 * (3.0 / 4.0 * kernel_2)

        self.__mem = []
        for i in range(0, self.__max_particle):
            p = Particle()
            p.identifier = i
            self.__mem.append(p)

        self.__cell = []
        for i in range(0, self.__total_cell):
            self.__cell.append(-1)

        self.__sys_running = 0

    def setup_system(self):
        print 'setup_system'
        velocity = [0.0, 0.0, 0.0]

        x = self.__world_size[0] * 0.2
        while x < self.__world_size[0] * 0.8:

            y = self.__world_size[1] * 0.5
            while y < self.__world_size[1] * 0.9:

                z = self.__world_size[2] * 0.2
                while z < self.__world_size[2] * 0.8:
                    self.__add_particle([x, y, z], velocity)
                    if self.__num_particle == self.__max_particle:
                        return

                    z += self.__kernel * 0.5
                y += self.__kernel * 0.5
            x += self.__kernel * 0.5

        return self.__num_particle

    def set_running_flag(self, running):
        self.__sys_running = running

    def update_particle(self):
        print 'update_particle'
        if self.__sys_running == 0:
            return

        self.__build_table()
        self.__compute_density_pressure()
        self.__compute_force_adv()
        self.__compute_advection()

    def get_mem(self):
        return self.__mem

    def __build_table(self):
        print 'build_table'
        for i in range(0, self.__total_cell):
            self.__cell[i] = -1

        for i in range(0, self.__num_particle):
            p = self.__mem[i]
            hash_value = self.__calc_cell_hash(self.__calc_cell_position(p.position))
            p.next_particle = self.__cell[hash_value]
            self.__cell[hash_value] = p.identifier

    def __compute_density_pressure(self):
        print 'compute_density_pressure'
        for i in range(0, self.__num_particle):
            p = self.__mem[i]
            cell_position = self.__calc_cell_position(p.position)

            p.density = 0.0
            p.pressure = 0.0

            for x in range(-1, 2):
                for y in range(-1, 2):
                    for z in range(-1, 2):
                        near_cell_position = [cell_position[0] + x, cell_position[1] + y, cell_position[2] + z]
                        hash_value = self.__calc_cell_hash(near_cell_position)

                        if hash_value == 0xffffffff:
                            continue

                        # 扫描kernel_2范围内的粒子
                        index = self.__cell[hash_value]
                        while index != -1:
                            np = self.__mem[index]
                            rel_position = [np.position[0] - p.position[0], np.position[1] - p.position[1],
                                            np.position[2] - p.position[2]]
                            r2 = (rel_position[0] * rel_position[0] + rel_position[1] * rel_position[1] +
                                  rel_position[2] * rel_position[2])

                            kernel_2 = self.__kernel * self.__kernel
                            if r2 < INF or r2 >= kernel_2:
                                index = np.next_index
                                continue

                            p.density = p.density + self.__mass * self.__poly6_value * math.pow(kernel_2 - r2, 3)
                            index = np.next_index

            # 密度及压力
            p.density = p.density + self.__self_density
            p.pressure = (math.pow(p.density / self.__rest_density, 7) - 1) * self.__gas_constant

    def __compute_force_adv(self):
        print 'compute_force_adv'
        for i in range(0, self.__num_particle):
            p = self.__mem[i]
            cell_position = self.__calc_cell_position(p.position)

            p.acceleration[0] = 0.0
            p.acceleration[1] = 0.0
            p.acceleration[2] = 0.0
            grad_color = [0.0, 0.0, 0.0]
            lplc_color = 0.0

            for x in range(-1, 2):
                for y in range(-1, 2):
                    for z in range(-1, 2):
                        near_position = [cell_position[0] + x, cell_position[1] + y, cell_position[2] + z]
                        hash_value = self.__calc_cell_hash(near_position)

                        if hash_value == 0xffffffff:
                            continue

                        index = self.__cell[hash_value]
                        while index != -1:
                            np = self.__mem[index]
                            rel_position = [p.position[0] - np.position[0], p.position[1] - np.position[1],
                                            p.position[2] - np.position[2]]
                            r2 = (rel_position[0] * rel_position[0] + rel_position[1] * rel_position[1] +
                                  rel_position[2] * rel_position[2])

                            kernel_2 = self.__kernel * self.__kernel
                            if kernel_2 > r2 > INF:
                                r = math.sqrt(r2)
                                v = self.__mass / np.density / 2.0
                                kernel_r = self.__kernel - r

                                # 计算压力项
                                pres_kernel = self.__spiky_value * kernel_r * kernel_r
                                temp_force = v * (p.pressure + np.pressure) * pres_kernel
                                p.acceleration[0] = p.acceleration[0] - rel_position[0] * temp_force / r
                                p.acceleration[1] = p.acceleration[1] - rel_position[1] * temp_force / r
                                p.acceleration[2] = p.acceleration[2] - rel_position[2] * temp_force / r

                                rel_vel = [np.ev[0] - p.ev[0], np.ev[1] - p.ev[1], np.ev[2] - p.ev[2]]
                                # 计算粘滞项
                                visc_kernel = self.__viscosity_value * (self.__kernel - r)
                                temp_force = v * self.__viscosity * visc_kernel
                                p.acceleration[0] = p.acceleration[0] + rel_vel[0] * temp_force
                                p.acceleration[1] = p.acceleration[1] + rel_vel[1] * temp_force
                                p.acceleration[2] = p.acceleration[2] + rel_vel[2] * temp_force

                                temp = -self.__grad_poly6 * v * pow(kernel_2 - r2, 2)
                                grad_color[0] += temp * rel_position[0]
                                grad_color[1] += temp * rel_position[1]
                                grad_color[2] += temp * rel_position[2]
                                lplc_color += self.__lplc_poly6 * v * (kernel_2 - r2) * (r2 - 3.0 / 4 * (kernel_2 - r2))

                            index = np.next_index

            lplc_color += self.__self_lplc_color / p.density
            p.surf_norm = math.sqrt(grad_color[0] * grad_color[0] + grad_color[1] * grad_color[1]
                                    + grad_color[2] * grad_color[2])

            if p.surf_norm > self.__surf_norm:
                p.acceleration[0] += self.__surf_coe * lplc_color * grad_color[0] / p.surf_norm
                p.acceleration[1] += self.__surf_coe * lplc_color * grad_color[1] / p.surf_norm
                p.acceleration[2] += self.__surf_coe * lplc_color * grad_color[2] / p.surf_norm

    def __compute_advection(self):
        print 'compute_advection'
        for i in range(0, self.__num_particle):
            p = self.__mem[i]
            # 速度及位置更新
            p.velocity[0] = (p.velocity[0] + p.acceleration[0] * self.__time_step / p.density +
                             self.__gravity[0] * self.__time_step)
            p.velocity[1] = (p.velocity[1] + p.acceleration[1] * self.__time_step / p.density +
                             self.__gravity[1] * self.__time_step)
            p.velocity[2] = (p.velocity[2] + p.acceleration[2] * self.__time_step / p.density +
                             self.__gravity[2] * self.__time_step)

            p.position[0] = p.position[0] + p.velocity[0] * self.__time_step
            p.position[1] = p.position[1] + p.velocity[1] * self.__time_step
            p.position[2] = p.position[2] + p.velocity[2] * self.__time_step

            # 墙壁
            if p.position[0] >= self.__world_size[0] - BOUNDARY:
                p.velocity[0] = p.velocity[0] * self.__wall_damping
                p.position[0] = self.__world_size[0] - BOUNDARY

            if p.position[0] < 0.0:
                p.velocity[0] = p.velocity[0] * self.__wall_damping
                p.position[0] = 0.0

            if p.position[1] >= self.__world_size[1] - BOUNDARY:
                p.velocity[1] = p.velocity[1] * self.__wall_damping
                p.position[1] = self.__world_size[1] - BOUNDARY

            if p.position[1] < 0.0:
                p.velocity[1] = p.velocity[1] * self.__wall_damping
                p.position[1] = 0.0

            if p.position[2] >= self.__world_size[2] - BOUNDARY:
                p.velocity[2] = p.velocity[2] * self.__wall_damping
                p.position[2] = self.__world_size[2] - BOUNDARY

            if p.position[2] < 0.0:
                p.velocity[2] = p.velocity[2] * self.__wall_damping
                p.position[2] = 0.0

            p.ev[0] = (p.ev[0] + p.velocity[0]) / 2.0
            p.ev[1] = (p.ev[1] + p.velocity[1]) / 2.0
            p.ev[2] = (p.ev[2] + p.velocity[2]) / 2.0

    def __add_particle(self, position, velocity):
        p = self.__mem[self.__num_particle]
        p.position[0] = position[0]
        p.position[1] = position[1]
        p.position[2] = position[2]

        p.velocity[0] = velocity[0]
        p.velocity[1] = velocity[1]
        p.velocity[2] = velocity[2]

        p.acceleration[0] = 0.0
        p.acceleration[1] = 0.0
        p.acceleration[2] = 0.0

        p.ev[0] = 0.0
        p.ev[1] = 0.0
        p.ev[2] = 0.0

        p.density = self.__rest_density
        p.pressure = 0.0

        self.__num_particle += 1

    def __calc_cell_position(self, position):
        cell_position = [0, 0, 0]
        cell_position[0] = int(math.floor(position[0] / self.__cell_size))
        cell_position[1] = int(math.floor(position[1] / self.__cell_size))
        cell_position[2] = int(math.floor(position[2] / self.__cell_size))
        return cell_position

    def __calc_cell_hash(self, cell_position):
        if cell_position[0] < 0 or cell_position[0] >= self.__grid_size[0] \
                or cell_position[1] < 0 or cell_position[1] >= self.__grid_size[1] \
                or cell_position[2] < 0 or cell_position[2] >= self.__grid_size[2]:
            return 0xffffffff

        pos = [0, 0, 0]
        pos[0] = cell_position[0] & (self.__grid_size[0] - 1)
        pos[1] = cell_position[1] & (self.__grid_size[1] - 1)
        pos[2] = cell_position[2] & (self.__grid_size[2] - 1)

        return (pos[2] * self.__grid_size[1] * self.__grid_size[0] +
                pos[1] * self.__grid_size[0] + pos[0])
