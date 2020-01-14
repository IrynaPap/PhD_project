from tkinter import messagebox
from calculations_in_one_point import calculations_velocities

profile_data = {}
profile_data_ex = {}


def get_data_from_file(data):
    data.pop(0)
    profile_data.update({'data': data})


def calculate_profile_depth_data():
    file = open("71918rightConc.csv", 'a')
    file.write("Q" + ' ' + "An" + ' ' + 'Bi' + ' ' + 'CPx' + ' ' + 'Ol' + ' ' + 'Amf' + ' ' + 'OPx' + ' '
               + 'Al' + ' ' + 'Ort' + ' ' + ' ' + 'h' + ' ' + 'x' + ' ' + 'y' + ' '
               + 'Vpinput' + ' ' + 'Vpoutput' + ' ' + 'Vs' + ' ' + 'Ro' + ' '
               + 'Sio2' + ' ' + "tiO2"  + ' ' + "al2O3" + ' ' + "fe2O3"  + ' ' + "feO"  + ' ' + "mnO" + ' ' + "mgO"  + ' '
               + "caO"  + ' ' + "na2O"  + ' ' + "k2O" + ' ' + "h2O" + '\n')
    file.close()
    pressure_temperature = {
        0: [0.0001, 12, 5.8],
        1: [0.0259, 19.83, 5.8],
        3: [0.0259, 19.83, 5.8],
        2: [0.0518, 28.65, 5.9],
        4: [0.1036, 53.2, 5.9],
        5: [0.1296, 70.9, 6],
        6: [0.1555, 93.1, 6],
        7: [0.1555, 93.1, 6],
        8: [0.2076, 146, 6.1],
        10: [0.2599, 200, 6.2],
        12: [0.3124, 245.04, 6.2],
        14: [0.3654, 281.764, 6.2],
        15: [0.392, 298, 6.3],
        16: [0.4186, 313.28, 6.3],
        18: [0.4886, 323.28, 6.3],
        20: [0.5267, 365, 6.4],
        22: [0.5867, 385, 6.5],
        24: [0.6267, 395, 6.5],
        25: [0.6659, 408, 6.6],
        30: [0.809, 432, 6.8]
    }
    # file = open("testfileProfileNEWWWWWWWWW.csv", 'a')
    #file.write('efVp, efVs, ro' + '\n')

    correct_real = {}
    if len(profile_data) > 0:
        for point in profile_data['data']:
            h = int(point[3])
            pressure = pressure_temperature[h][0]
            temperature = pressure_temperature[h][1]
            velocity = float(point[4])
            x = float(point[0])
            y = float(point[1])
            correct_realisations_sio2 = calculations_velocities(pressure, temperature, velocity, h, x, y)
            if len(correct_realisations_sio2) != 0:
                correct_real.update({h: correct_realisations_sio2})
            # print(h)
            #file.write(str(h) + ' ' + str(correct_realisations_sio2) + '\n')
        return correct_real
    else:
        messagebox.showwarning("Warning", "Download profile data")


def send_profile_data():
    return profile_data