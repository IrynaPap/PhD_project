# from form import data_for_downloading_profile_import
import random
import math
import cmath

# TO DO вынести минерал и температура давление в отдельный класс
fn = 0
filename = "all_realisations" + str(fn)
all_realizations_file = open(str(filename) + ".csv", "w")

file = open("clc_profile_8.csv", 'a')

# mineral_gl = {
#     "Q": {"C1": 0.01, "F1": 0.999, "k": 39.5, "G": 48.8, "dk/dP": 4.3, "dk/dT": -0.003, "dG/dP": -0.45, "dG/dT": -0.008,
#           "alfa": 5, "ro": 2.646, "dLr/dP": -0.8, "dLr/dT": 0.007, "siO2": 100},
#     "An": {"C1": 0.99, "F1": 0.999, "k": 88.6, "G": 41, "dk/dP": 4.3, "dk/dT": -0.015, "dG/dP": 3.9, "dG/dT": -0.0056,
#            "alfa": 8.6, "ro": 2.77, "dLr/dP": -0.57, "dLr/dT": 0.008, "siO2": 44.2714},
#     "Bi": {"C1": 0, "F1": 0.999, "k": 58.5, "G": 24.1, "dk/dP": 5.43, "dk/dT": -0.003, "dG/dP": 3.82, "dG/dT": -0.0016,
#            "alfa": 9.59, "ro": 2.905, "dLr/dP": -0.56, "dLr/dT": 0.0048, "siO2": 36.2038},
#     "CPx": {"C1": 0, "F1": 0.999, "k": 112.9, "G": 67.1, "dk/dP": 4, "dk/dT": -0.003, "dG/dP": 2, "dG/dT": -0.0019,
#             "alfa": 9.7, "ro": 3.28, "dLr/dP": -0.7, "dLr/dT": 0.0008, "siO2": 50.63},
#     "Ol": {"C1": 0, "F1": 0.999, "k": 131.23, "G": 76.49, "dk/dP": 5.38, "dk/dT": -0.018, "dG/dP": 1.76,
#            "dG/dT": -0.0135, "alfa": 9.65, "ro": 3.27, "dLr/dP": -1.416, "dLr/dT": 0.0035, "siO2": 39.4029},
#     "Amf": {"C1": 0, "F1": 0.999, "k": 92.4, "G": 50.4, "dk/dP": 2.41, "dk/dT": -0.0063, "dG/dP": 3.7, "dG/dT": -0.0111,
#             "alfa": 11.2, "ro": 3.041, "dLr/dP": -0.56, "dLr/dT": 0.0033, "siO2": 45.531},
#     "OPx": {"C1": 0, "F1": 0.999, "k": 107.5, "G": 72.6, "dk/dP": 9.61, "dk/dT": -0.037, "dG/dP": 2.18,
#             "dG/dT": -0.0166, "alfa": 15.7, "ro": 3.584, "dLr/dP": -3.023, "dLr/dT": 0.0024, "siO2": 48.9935},
#     "Al": {"C1": 0, "F1": 0.999, "k": 57.9, "G": 32.3, "dk/dP": 3.91, "dk/dT": -0.008, "dG/dP": 1.3, "dG/dT": -0.0051,
#            "alfa": 11.4, "ro": 2.552, "dLr/dP": -0.74, "dLr/dT": 0.0065, "siO2": 67.4715},
#     "Ort": {"C1": 0, "F1": 0.999, "k": 59.8, "G": 30.1, "dk/dP": 4.22, "dk/dT": -0.0075, "dG/dP": 2.3, "dG/dT": -0.0045,
#             "alfa": 8.22, "ro": 2.546, "dLr/dP": -0.49, "dLr/dT": 0.0054, "siO2": 67.4715}
# }

mineral_gl = {
    "Q": {"C1": 0.01, "F1": 0.999, "k": 39.5, "G": 48.8, "dk/dP": 4.3, "dk/dT": -0.003, "dG/dP": -0.45, "dG/dT": -0.008,
          "alfa": 5, "ro": 2.646, "dLr/dP": -0.8, "dLr/dT": 0.007, "siO2": 100,
          "tiO2": 0, "al2O3": 0, "fe2O3": 0, "feO": 0, "mnO": 0, "mgO": 0,
          "caO": 0, "na2O": 0, "k2O": 0, "h2O": 0},

    "An": {"C1": 0.99, "F1": 0.999, "k": 88.6, "G": 41, "dk/dP": 4.3, "dk/dT": -0.015, "dG/dP": 3.9, "dG/dT": -0.0056,
           "alfa": 8.6, "ro": 2.77, "dLr/dP": -0.57, "dLr/dT": 0.008, "siO2": 44.2714,
           "tiO2": 0.025, "al2O3": 35.4186, "fe2O3": 0.4283, "feO": 0.02666, "mnO": 0, "mgO": 0.115,
           "caO": 18.6629, "na2O": 0.5942, "k2O": 0.12, "h2O": 0.4971},

    "Bi": {"C1": 0, "F1": 0.999, "k": 58.5, "G": 24.1, "dk/dP": 5.43, "dk/dT": -0.003, "dG/dP": 3.82, "dG/dT": -0.0016,
           "alfa": 9.59, "ro": 2.905, "dLr/dP": -0.56, "dLr/dT": 0.0048, "siO2": 36.2038,
           "tiO2": 2.7652, "al2O3": 16.0504, "fe2O3": 4.4599, "feO": 16.0518, "mnO":0.2626, "mgO": 10.3228,
           "caO":1.0633, "na2O": 0.4175, "k2O": 8.2127, "h2O": 2.077},

    "CPx": {"C1": 0, "F1": 0.999, "k": 112.9, "G": 67.1, "dk/dP": 4, "dk/dT": -0.003, "dG/dP": 2, "dG/dT": -0.0019,
            "alfa": 9.7, "ro": 3.28, "dLr/dP": -0.7, "dLr/dT": 0.0008, "siO2": 50.63,
            "tiO2": 0.4, "al2O3": 2.08, "fe2O3": 2.25, "feO": 9.19, "mnO": 0.34, "mgO": 12.16,
            "caO": 21.11, "na2O": 0.64, "k2O": 0.08, "h2O": 0.12},

    "Ol": {"C1": 0, "F1": 0.999, "k": 131.23, "G": 76.49, "dk/dP": 5.38, "dk/dT": -0.018, "dG/dP": 1.76,
           "dG/dT": -0.0135, "alfa": 9.65, "ro": 3.27, "dLr/dP": -1.416, "dLr/dT": 0.0035, "siO2": 39.4029,
           "tiO2": 0.0645, "al2O3": 0.3043, "fe2O3": 1.0111, "feO": 15.0999, "mnO": 0.4703, "mgO": 43.7496,
           "caO": 0.3795, "na2O": 0.10511, "k2O": 0.0438, "h2O": 0.182},

    "Amf": {"C1": 0, "F1": 0.999, "k": 92.4, "G": 50.4, "dk/dP": 2.41, "dk/dT": -0.0063, "dG/dP": 3.7, "dG/dT": -0.0111,
            "alfa": 11.2, "ro": 3.041, "dLr/dP": -0.56, "dLr/dT": 0.0033, "siO2": 45.531,
            "tiO2": 1.1909, "al2O3": 10.3716, "fe2O3": 4.0457, "feO": 11.6334, "mnO": 0.2867, "mgO": 9.461,
            "caO": 10.6416, "na2O": 1.534, "k2O": 0.6645, "h2O": 2.8617},

    "OPx": {"C1": 0, "F1": 0.999, "k": 107.5, "G": 72.6, "dk/dP": 9.61, "dk/dT": -0.037, "dG/dP": 2.18,
            "dG/dT": -0.0166, "alfa": 15.7, "ro": 3.584, "dLr/dP": -3.023, "dLr/dT": 0.0024, "siO2": 48.9935,
            "tiO2": 0.2421, "al2O3": 1.5561, "fe2O3": 3.2499, "feO": 29.2095, "mnO": 0.5038, "mgO": 15.4605,
            "caO": 0.9638, "na2O": 0.1812, "k2O": 0.1219, "h2O": 0.31},

    "Al": {"C1": 0, "F1": 0.999, "k": 57.9, "G": 32.3, "dk/dP": 3.91, "dk/dT": -0.008, "dG/dP": 1.3, "dG/dT": -0.0051,
           "alfa": 11.4, "ro": 2.552, "dLr/dP": -0.74, "dLr/dT": 0.0065, "siO2": 67.4715,
           "tiO2": 0.005, "al2O3": 20.0426, "fe2O3": 0.2358, "feO": 0, "mnO": 0.045, "mgO": 0.17,
           "caO": 0.46678, "na2O": 11.0723, "k2O": 0.4298, "h2O": 0.3712},

    "Ort": {"C1": 0, "F1": 0.999, "k": 59.8, "G": 30.1, "dk/dP": 4.22, "dk/dT": -0.0075, "dG/dP": 2.3, "dG/dT": -0.0045,
            "alfa": 8.22, "ro": 2.546, "dLr/dP": -0.49, "dLr/dT": 0.0054, "siO2": 64.1052,
            "tiO2": 0.0108, "al2O3": 19.3359, "fe2O3": 0.2426, "feO": 0.0875, "mnO": 0.0033, "mgO": 0.072,
            "caO": 0.5278, "na2O": 3.0103, "k2O": 12.2991, "h2O": 0.3331},
    #granat
}

global count_matrix
count_matrix = {"realisations_count": 2,
                "matrix_mineral": "Al"}


def get_data(data):
    count_matrix["realisations_count"] = int(data[1])
    count_matrix["matrix_mineral"] = data[0]




count = count_matrix["realisations_count"]



def calculations_velocities(pressure, temperature, velocity, h, x, y):
    correct_realisations_sio2 = []
    single_realisation = {}
    mineral = mineral_gl

    def shuffle_min_components():
        shuffle_min_comp = []
        for component in mineral:
            shuffle_min_comp.append(component)
        random.shuffle(shuffle_min_comp)
        return shuffle_min_comp

    # Берем из ассоциативного массива miner элементы(минералы) и случайным образом перемешиваем

    def rand_components_percent(shuffle_min_comp):
        max_comp = 1
        sum_comp = 0
        i = 0
        massive_length = len(shuffle_min_comp)
        for component in shuffle_min_comp:
            if i < massive_length - 1:
                rand_comp_c = random.uniform(0, max_comp)
                mineral[component]["C1"] = rand_comp_c
            else:
                mineral[component]["C1"] = max_comp
            max_comp = max_comp - rand_comp_c
            sum_comp = rand_comp_c + sum_comp
            i = i + 1
        return mineral

    correct_realisations = []

    # Случайным образом генерируем концентрации/ На вход функции передаем перемашаный массив компонентов и записываем концентрации сразу в массив miner
    # print(count_matrix["realisations_count"],  count_matrix["matrix_mineral"])
    # t = pressure_temperature[1][0]

    def calculate_minerals_properties(mineral):

        for min_comp in mineral:
            elast2 = mineral[min_comp]["k"] + mineral[min_comp]["dk/dP"] * pressure + mineral[min_comp][
                "dk/dT"] * temperature  # {kk1}
            mineral[min_comp].update({"Elast2": elast2})
            elast3 = mineral[min_comp]["G"] + mineral[min_comp]["dG/dP"] * pressure + mineral[min_comp][
                "dG/dT"] * temperature  # {mu1}
            mineral[min_comp].update({"Elast3": elast3})
            elast4 = (mineral[min_comp]["alfa"] + mineral[min_comp]["dLr/dP"] * pressure + mineral[min_comp][
                "dLr/dT"] * temperature) * 1E-6  # {Lr} //{Add 1E-6 [1/k]}
            mineral[min_comp].update({"Elast4": elast4})
            ro_pt = mineral[min_comp]["ro"] * (1 + pressure / elast2 - 3 * elast4 * temperature)
            mineral[min_comp].update({"roPT": ro_pt})
            elast5 = elast4 * 3 * elast2 * 1E+9
            mineral[min_comp].update({"Elast5": elast5})
        # что то может пойти не так
        for min_comp in mineral:
            matrix_min_elast5 = mineral[count_matrix["matrix_mineral"]]["Elast5"]
            elast6 = mineral[min_comp]["Elast5"] - matrix_min_elast5
            # print(min_comp, mineral[min_comp]["Elast5"], matrix_min_elast5)
            mineral[min_comp].update({"Elast6": elast6})
        return mineral

    def calculate_density(mineral):
        ro = 0
        for key in mineral:
            ro = ro + mineral[key]["roPT"] * mineral[key]["C1"]
            single_realisation.update({"ro": ro})
            single_realisation.update({"mineral_exemplar": mineral})
        # print(ro)
        return single_realisation

    def calculate_cm_matrix():
        cm = mineral[count_matrix["matrix_mineral"]]["C1"]
        return cm

    def calculate_current_elast_modules(realisation, cm):
        m0 = 0
        k0 = 0
        j1 = {}
        r = {}
        c = {};
        m1 = {};
        k1 = {}
        j = {};
        j1 = {};
        j2 = {};
        j3 = {};
        r3 = {};
        kh = {}
        kg = {};
        kgt = {};
        mg = {};
        mgt = {};
        mf = {};
        kf = {};
        mc = {};
        mct = {};
        kc = {};
        kct = {};
        tz = {};
        lz = {};
        kz = {};
        nz = {}
        vc = {};
        wc = {};
        tc = {};
        lc = {};
        lct = {};
        nc = {};
        cc = {};
        lzt = {};
        vz = {};
        wz = {};
        mz = {};
        md = {};
        kd = {};
        mh = {};
        l1 = {};
        nh = {};
        lh = {};
        fh = {};
        u1 = {};
        nu2 = {};
        nu3 = {};
        tg = {};
        lg = {};
        ng = {};
        vg = {};
        wg = {}
        mineral_exemplar = realisation["mineral_exemplar"]
        #  density = realisation["ro"]
        for min_comp in mineral_exemplar:
            m1[min_comp] = 2 * mineral_exemplar[min_comp]["Elast3"]
            k1[min_comp] = 3 * mineral_exemplar[min_comp]["Elast2"]
            m0 = m0 + mineral_exemplar[min_comp]["C1"] * m1[min_comp]
            k0 = k0 + mineral_exemplar[min_comp]["C1"] * k1[min_comp]

        t0 = (2 * k0 + m0) / 3
        n0 = (k0 + 2 * m0) / 3
        g0 = (-1) / (2 * m0)
        g1 = (-1) / (2 * n0)
        # k = -1
        for mi in mineral_exemplar:
            # k = k + 1
            min_mi = mineral[mi]
            r = (mineral[mi]["F1"] ** 2) - 1
            c = math.sqrt(math.fabs(r))
            if mineral[mi]["F1"] > 1:
                j = math.log(mineral[mi]["F1"] + c)
            elif mineral[mi]["F1"] < 1:
                j = math.atan(c / mineral[mi]["F1"])
            else:
                j = 1.5

            j1[mi] = (j * mineral[mi]["F1"] / c - 1) / r
            j2[mi] = 1 - j1[mi]
            j3[mi] = ((1 + 2 * math.pow(mineral[mi]["F1"], 2)) * j1[mi] - 1) / (2 * r)
            r3[mi] = j3[mi] * t0 / m0
            tg[mi] = g1 * (j2[mi] + r3[mi])
            lg[mi] = (-g1) * r3[mi]
            ng[mi] = 2 * g1 * (j1[mi] + r3[mi])
            vg[mi] = g0 * j2[mi] + (tg[mi] / 2)
            wg[mi] = g0 * (1 + j1[mi]) + 2 * lg[mi]
            kg[mi] = lg[mi] + tg[mi]
            kgt[mi] = ng[mi] + 2 * lg[mi]
            mg[mi] = tg[mi] - 2 * lg[mi]
            mgt[mi] = ng[mi] - lg[mi]
            mf[mi] = m1[mi] - m0
            kf[mi] = k1[mi] - k0
            mc[mi] = 1 - mf[mi] * mg[mi]
            mct[mi] = 1 - mf[mi] * mgt[mi]
            kc[mi] = 1 - kf[mi] * kg[mi]
            kct[mi] = 1 - kf[mi] * kgt[mi]
            vc[mi] = 1 - mf[mi] * vg[mi]
            wc[mi] = 1 - mf[mi] * wg[mi]
            tc[mi] = (2 * kc[mi] + mc[mi]) / 3
            lc[mi] = (kc[mi] - mc[mi]) / 3
            lct[mi] = (kct[mi] - mct[mi]) / 3
            nc[mi] = (kct[mi] + 2 * mct[mi]) / 3
            cc[mi] = nc[mi] * tc[mi] - 2 * lc[mi] * lct[mi]
            tz[mi] = (nc[mi] * tg[mi] - 2 * lc[mi] * lg[mi]) / cc[mi]
            lz[mi] = (nc[mi] * lg[mi] - lc[mi] * ng[mi]) / cc[mi]
            lzt[mi] = (tc[mi] * lg[mi] - lct[mi] * tg[mi]) / cc[mi]
            nz[mi] = (tc[mi] * ng[mi] - 2 * lct[mi] * lg[mi]) / cc[mi]
            vz[mi] = vg[mi] / vc[mi]
            wz[mi] = wg[mi] / wc[mi]
            mz[mi] = 1 / 15 * (tz[mi] + 2 * nz[mi] - 2 * lz[mi] - 2 * lzt[mi] + 6 * vz[mi] + 6 * wz[mi])
            kz[mi] = 1 / 3 * (2 * (tz[mi] + lz[mi] + lzt[mi]) + nz[mi])
            mhm0 = 0
            khm0 = 0

        for mi2 in mineral:
            km = k1[count_matrix["matrix_mineral"]]
            mm = m1[count_matrix["matrix_mineral"]]
            md[mi2] = m1[mi] - mm
            kd[mi2] = k1[mi] - km
            mhm0 = mhm0 + min_mi["C1"] * md[mi2] * mz[mi2]
            khm0 = khm0 + min_mi["C1"] * kd[mi2] * kz[mi2]
        mhm = 1 / (1 + mhm0)
        khm = 1 / (1 + khm0)
        for mi1 in mineral:
            mh[mi1] = mhm * (1 + md[mi1] * mz[mi1])
            kh[mi1] = khm * (1 + kd[mi1] * kz[mi1])
        kk0 = 0
        mm0 = 0
        for mi3 in mineral:
            min_mi3 = mineral[mi3]
            kk0 = kk0 + min_mi3["C1"] * kh[mi3] * k1[mi3]
            mm0 = mm0 + min_mi3["C1"] * mh[mi3] * m1[mi3]

        kv = (kk0 + cm * km * khm) / 3
        mu = (mm0 + cm * mm * mhm) / 2
        # print("kv,mu", h, kv, mu)
        # где то прошляпились десятые
        single_realisation.update({"Elast0": kv})
        single_realisation.update({"Elast1": mu})
        # вызываем функцию calc_effective_linear и передаем kz, khm, km, k1, kh

        calc_effective_linear(kz, khm, km, k1, kh)

    # Rock = single_realisation
    # TO DO
    def calc_effective_linear(kz, khm, km, k1, kh):
        alfa_r = []
        alfa_rm = 0
        ss_beta = 0
        ef_beta = 0
        ii = -1
        r_beta = 0

        mineral_keys = list(mineral)
        cm = calculate_cm_matrix()

        for mi in mineral:
            min_mi = mineral[mi]
        if mi != count_matrix["matrix_mineral"]:
            ii = ii + 1
            r_beta = r_beta + mineral[mi]["C1"] * kz[mineral_keys[ii]] * min_mi["Elast6"] * 1E-9

        ii = -1
        for mi in mineral:
            min_mi = mineral[mi]
            if mi != count_matrix["matrix_mineral"]:
                ii = ii + 1
                kh1 = kh[mineral_keys[ii]]
                kz1 = kz[mineral_keys[ii]]
                alfa_r_i = kh1 * r_beta - kz1 * min_mi["Elast6"] * 1E-9
                alfa_r.append(alfa_r_i)
            alfa_rm = khm * r_beta
            ss_beta = 0

        ii = -1
        for mi in mineral:
            min_mi = mineral[mi]
            if mi != count_matrix["matrix_mineral"]:
                ii = ii + 1
                ss_beta = ss_beta + mineral[mi]["C1"] * min_mi["Elast5"] - alfa_r[ii] * k1[mineral_keys[ii]]
        ef_beta = ss_beta + cm * (mineral[count_matrix["matrix_mineral"]]["Elast5"] - alfa_rm * km)
        fef_lr = ef_beta / (3 * single_realisation["Elast0"] * 1E+9)
        ef_lr = fef_lr
        # вызываем функцию calc_effective_kg и передаем давление температуру и значение ef_lr
        calc_effective_kg(ef_lr)

    # calcEffectiveKG(0.6659, 408)

    def calc_effective_kg(ef_lr):
        velocities = {}
        solutions_vector = []
        nu = (3 * single_realisation["Elast0"] - 2 * single_realisation["Elast1"]) / (
                    6 * single_realisation["Elast0"] + 2 * single_realisation["Elast1"])
        rp = single_realisation["ro"] * (1 + pressure / single_realisation["Elast0"] - 3 * ef_lr * temperature * 1E-6)
        # print("rp", rp, single_realisation["ro"])
        # ef_lr1 зачем нужно
        ef_lr1 = (1 - math.sqrt(nu)) / (1 - 2 * nu)
        ro = rp
        ef_k = single_realisation["Elast0"]
        ef_g = single_realisation["Elast1"]
        ef_vp = math.sqrt((ef_k + 4 / 3 * ef_g) / ro)
        ef = ((ef_k + 4 / 3 * ef_g) / ro)
        ef_vs = math.sqrt(ef_g / ro)
        velocities.update({"h": h})
        velocities.update({"efVp": ef_vp})
        velocities.update({"efVs": ef_vs})
        velocities.update({"ro": ro})
        velocities.update({"mineral": mineral})
        solutions_vector.append(velocities)
        # m += m непонятно что
        # print(h, "efvp", ef_vp, ro) #, pressure, temperature, ef_k, ef_g, "ro", ro)
        # print(solutions_vector)
        # print(solutions_vector[0]['h'], solutions_vector[0]['efVp'], solutions_vector[0]['efVs'], solutions_vector[0]['ro'],
        #     solutions_vector[0]['mineral']['Q']['C1'], solutions_vector[0]['mineral']['An']['C1'],
        #    solutions_vector[0]['mineral']['Bi']['C1'], solutions_vector[0]['mineral']['CPx']['C1'],
        #   solutions_vector[0]['mineral']['Ol']['C1'], solutions_vector[0]['mineral']['OPx']['C1'],
        #  solutions_vector[0]['mineral']['Amf']['C1'], solutions_vector[0]['mineral']['Al']['C1'])

        select_right_realisations(solutions_vector)

    def sio_two_concentration(correct_realisations, oxyd):

        # print(correct_realisations)
        w = 0
        m = 0
        ro_min = 0
        c_min = 0
        sio2_rock = 0
        oxyd_sum = 0
        mineral_keys = list(mineral)
        m_sum = 0
        for l in range(0, len(mineral_keys)): # по хім композіціі
            # ro_min = mineral[mineral_keys[l]]["ro"]

            oxyd_conc = mineral[mineral_keys[l]][oxyd]
            c_min = mineral[mineral_keys[l]]["C1"]
            m = oxyd_conc * c_min
            # print(oxyd_conc, c_min, m, oxyd_sum)
            # m_rock.append(m)
            # print(oxyd_conc, c_min)
            oxyd_sum += m
            file = open("right14Conc.csv", 'a')
            file.write( str(c_min) + ' ')
            #имена миенралов, оксидов, скорости, плостности, x, y, h
        file.write('\n')
        file.close()
        print(h, oxyd_sum)
            # file = open(str(oxyd + ".csv"), 'a')
            # file.write(str(correct_realisations['h'])+ ' '+str(correct_realisations['mineral'])+' '+str(sio2_rock)+"\n")
            # file.close()


    # Создаем вектор реализаций
    # correct_realisations = []
    def select_right_realisations(vector_of_realisations):
        velocity_etalon = velocity
        # print(vector_of_realisations[0]['h'], vector_of_realisations[0]['efVp'], vector_of_realisations[0]['ro'], vector_of_realisations[0]['mineral']['Q']['C1'])
        # all_realizations_file.write("A"+" "+str(vector_of_realisations[0]['h'])+" "+str(vector_of_realisations[0]['efVp'])+" "+str(vector_of_realisations[0]['efVs'])+" "+str(vector_of_realisations[0]['ro'])+" "+str(vector_of_realisations[0]['mineral']['Q']['C1'])+" "+str(vector_of_realisations[0]['mineral']['An']['C1'])+" "+str(vector_of_realisations[0]['mineral']['Bi']['C1'])+" "+str(vector_of_realisations[0]['mineral']['CPx']['C1'])+" "+str(vector_of_realisations[0]['mineral']['Ol']['C1'])+" "+str(vector_of_realisations[0]['mineral']['OPx']['C1'])+" "+str(vector_of_realisations[0]['mineral']['Amf']['C1'])+" "+str(vector_of_realisations[0]['mineral']['Al']['C1'])+'\n')
        correct_realisations = []

        vector_length = len(vector_of_realisations)
        # print(vector_of_realisations)
        for i in range(0, vector_length):
            solutions_vector_value = math.fabs(vector_of_realisations[i]["efVp"])
            value = (round(solutions_vector_value, 2) - velocity_etalon) / velocity_etalon
            value_round = round(value, 2)
            if math.fabs(value_round) <= 0.03:
                correct_realisations.append(vector_of_realisations[i])
                file = open("rightConc.csv", 'a')
                # file.write(str(velocity_etalon) + ' ' + str(vector_of_realisations[i]['C1']) + '\n')
                file.close()
                sio_two_concentration(correct_realisations[i], "siO2")
                # print(vector_of_realisations[i])
            else:
                m = "graite sche"

        # sio_two_concentration(correct_realisations, "na2O")
        # sio_two_concentration(correct_realisations, "tiO2")
        # sio_two_concentration(correct_realisations, "fe2O3")
        # sio_two_concentration(correct_realisations, "al2O3")
        # sio_two_concentration(correct_realisations, "feO")
        # sio_two_concentration(correct_realisations, "mnO")
        # sio_two_concentration(correct_realisations, "mgO")



    def solutions_vector_generation():
        efVP = 0
        efVs = 0
        ro = 0
        for i in range(0, count_matrix["realisations_count"]):
            shuffle_min_component = shuffle_min_components()
            rand_component_percent = rand_components_percent(shuffle_min_component)
            calculate_mineral_properties = calculate_minerals_properties(rand_component_percent)
            cm = calculate_cm_matrix()
            realisation_with_density = calculate_density(calculate_mineral_properties)
            calculate_current_elast_modules(realisation_with_density, cm)
            # print(realisation_with_densitysity)

        # #print(h, len(correct_realisations_sio2), correct_realisations_sio2)
        # for i in range(0, len(correct_realisations_sio2)):
        #     efVp = correct_realisations_sio2[i][0]['efVp']
        #     efVs = correct_realisations_sio2[i][0]['efVs']
        #     ro = correct_realisations_sio2[i][0]['ro']
        #     mineral = correct_realisations_sio2[i][0]['mineral']
        #     Q = str(mineral['Q']['C1'])
        #     Bi = str(mineral['Bi']['C1'])
        #     An = str(mineral['An']['C1'])
        #     Cpx = str(mineral['CPx']['C1'])
        #     Ol = str(mineral['Ol']['C1'])
        #     Amf = str(mineral['Amf']['C1'])
        #     Opx = str(mineral['OPx']['C1'])
        #     Al = str(mineral['Al']['C1'])
        #     # print(mineral)
        #     # print(h, i, Q, Bi)
        #     #print(h, len(correct_realisations_sio2), efVp, ro)
        #     #file = open("8ExportHXYEfVp.csv", 'a')
        #     #file.write(str(h) + ' ' + str(x) + ' ' + str(y) + ' ' + str(efVp) + '\n')
        #     #file.close()
        #     #file = open("8ExportMineralComp.csv", 'a')
        #     #file.write(str(h) + ' ' + str(x) + ' ' + str(y) + ' ' + Q + ' ' + Bi + ' ' + An + ' ' + Cpx + ' ' + Ol + ' ' + Amf + ' ' + Opx + ' ' + Al + '\n')
        #     #file.close()

        #print(h, correct_realisations_sio2)

    solutions_vector_generation()
    return correct_realisations_sio2
