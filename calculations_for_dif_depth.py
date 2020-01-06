from calculations_in_one_point import calculations_velocities


def calculate_on_the_depth():

    pressure_temperature = {
            0: [0.0001, 12, 5.8],
            5: [0.1296, 70.9, 6],
            10: [0.2599, 200, 6.15],
            20: [0.5267, 365, 6.4],
            30: [0.809, 432, 6.8]
        }
    correct_real = {}
    real_h = []
    i = 0
    for h in pressure_temperature:
        pressure = pressure_temperature[h][0]*1.1
        temperature = pressure_temperature[h][1]
        velocity = pressure_temperature[h][2]
        correct_realisations_sio2 = calculations_velocities(pressure, temperature, velocity, h)
        if len(correct_realisations_sio2) != 0:
            real_h.append(correct_realisations_sio2)
            #print(h, real_h)
            #print(correct_realisations_sio2)
            correct_real.update({h: correct_realisations_sio2})
           # correct_realisations_sio23.append(correct_realisations_sio2)
            #correct_realisations_sio23.append(h)
            i=i+1
    #print(correct_realisations_sio2)
    return correct_real



