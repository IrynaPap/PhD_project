#text file parsing
import re
#import pandas as pd
from calculation_for_profile_depth import get_data_from_file

data_from_file = {}
col_data = []


def parse_text_data(file):
        col_num = 2

        delimiter = "\t"
        delimiter2 = "\n"
        with open(file, 'r') as f:
            for line in f:
                line2 = line.replace("\n", "")
                line3 = line2.replace(",", ".")
                col_data.append(line3.split(delimiter))
        get_data_from_file(col_data)
        #data_from_file.update({"data":col_data})

    #print(data_from_file)


#def send_data_from_file():
#    return data_from_file