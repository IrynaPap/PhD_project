from calculation_for_profile_depth import calculate_profile_depth_data
from calculations_for_dif_depth import calculate_on_the_depth



def get_calculated_profile_data():
    data = calculate_profile_depth_data()
    # data1 = str(send_profile_data())
    file = open("testfileProfile.csv", "w")
    file.write('efVp, efVs, ro' + '\n')
    for d in data:
        #print(d)
        #print(data[d])
        for r in data[d]:
            h = str(r[0]['h'])
            efVp = str(r[0]['efVp'])
            efVs = str(r[0]['efVs'])
            ro = str(r[0]['ro'])
        file.write(str(d) + ' ' + h + ' ' + efVp + ' ' + efVs + ' ' + ro + '\n')
    file.close()
    #print(data)
    return data

def get_calculated_default_data():
    data = calculate_on_the_depth()
    conc = []
    file = open("testfile.csv", "w")
    file.write('efVp, efVs, ro, sio2, Q, Bi, An, Cpx, Ol, Amf, Opx, Al'+'\n')
    #print(data)
    for d in data:
     #   print(d)

        for r in data[d]:
            efVp = str(r[0]['efVp'])
            efVs = str(r[0]['efVs'])
            ro = str(r[0]['ro'])
            sio2 = str(r[0]['sio2'])
            mineral = r[0]['mineral']
            Q = str(mineral['Q']['C1'])
            Bi = str(mineral['Bi']['C1'])
            An = str(mineral['An']['C1'])
            Cpx = str(mineral['CPx']['C1'])
            Ol = str(mineral['Ol']['C1'])
            Amf = str(mineral['Amf']['C1'])
            Opx = str(mineral['OPx']['C1'])
            Al = str(mineral['Al']['C1'])
        file.write(str(d)+' '+efVp + ' ' + efVs + ' ' + ro + ' ' + sio2 + ' ' + Q + ' ' + Bi + ' ' + An + ' ' + Cpx + ' ' + Ol + ' ' + Amf + ' ' + Opx + ' ' + Al +'\n')
    file.close()
