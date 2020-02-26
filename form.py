from draw_plot import drawing_plots_by_call
from draw_plot import draw_profiles
from data_parsing import parse_text_data
from calculations_in_one_point import get_data
from calculations_for_dif_depth import calculate_on_the_depth
from calculation_for_profile_depth import calculate_profile_depth_data
from create_data_for_file import get_calculated_profile_data
from create_data_for_file import get_calculated_default_data
from calculations_in_one_point import *
from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfilename


root = Tk()
data_import_internal = {}


def create_form():
    #fn = 0
    def get_input_data():
        matr_and_count = []
        matrix_el = mtrix_element.get()
        count_iter = count.get()
        matr_and_count.append(matrix_el)
        matr_and_count.append(count_iter)
        # matr_and_count.append(8)
        get_data(matr_and_count)

    def draw_plot():
        get_input_data()
        calc_real = calculate_on_the_depth()
        drawing_plots_by_call(calc_real)
        filename = "all_realisations" + str(fn)
        all_realizations_file = open(str(filename) + ".csv", "w")


    def calculate_profile_data():
        get_input_data()
        calc_real = calculate_profile_depth_data()
        # draw_profiless(calc_real)

    def download_data():
        data = get_calculated_profile_data()
        file = open("testfileNew1.csv", "w")
        print('loading works')
        file.write(str(data))
        file.close()

    def download_default_data():
        get_calculated_default_data()

    def load_file():
        profile_data = askopenfilename(initialdir="/Users/irynapap/Documents/PhD/Исходники_код_Паскаль/CMC",
                               filetypes=(("Text File", "*.txt"), ("All Files", "*.*")),
                               title="Choose a file."
                               )
        try:
            parse_text_data(profile_data)
            load_completed = Label(text="File is loaded")
            load_completed.pack()
            load_completed.place(relx=.7, rely=.2)
            return profile_data
        except:
            print("No file exists")

    root.title("Calculate mineral composition")
    root.geometry("600x400+300+350")
    count = StringVar()

    count_of_iterations = Entry(textvariable=count)
    count_of_iterations.place(relx=.3, rely=.1, anchor="c")
    count_of_iterations.insert(0, "1")

    # profile_number = Entry(textvariable=count)
    # profile_number.place(relx=.3, rely=.4, anchor="c")
    # profile_number.insert(8, "8")


    OPTIONS = [
        "An",
        "Q",
        "Al",
        "Bi",
        "OPx",
        "Amf",
        "CPx"
    ]


    mtrix_element = StringVar(root)
    mtrix_element.set(OPTIONS[0])  # default value
    w = OptionMenu(root, mtrix_element, *OPTIONS)
    w.pack()

    w.place(relx=.3, rely=.3, anchor="c")

    btn = Button(text="Calculate default data", command=draw_plot)
    btn.place(relx=.5, rely=.7, anchor="c")

    btn_load_file = Button(text="Load profile data", command=load_file)
    btn_load_file.place(relx=.8, rely=.1, anchor="c")

    btn_download_file = Button(text="Download calculated profile", command=download_data)
    btn_download_file.place(relx=.8, rely=.3, anchor="c")

    btn_download_default_file = Button(text="Download calculated default", command=download_default_data)
    btn_download_default_file.place(relx=.8, rely=.4, anchor="c")

    btn_load_file = Button(text="Calculate profile data", command=calculate_profile_data)
    btn_load_file.place(relx=.8, rely=.2, anchor="c")

    root.mainloop()


create_form()
#Это то что вызывает рисовалку плотов
#drawing_plots_by_call()



