from tkinter import *
from tkinter import ttk

import numpy as np  

def	Start(*args):
    try:
        # Pull Inputs
        x_dim = float(X_input.get())
        y_dim = float(Y_input.get())
        z_dim = float(Z_input.get())

        k_const = float(k_input.get())

        input_west = str(bc_input_west.get())
        input_north = str(bc_input_north.get())
        input_east = str(bc_input_east.get())
        input_south = str(bc_input_south.get())

        value_west = float(bc_value_west.get())
        value_north = float(bc_value_north.get())
        value_east = float(bc_value_east.get())
        value_south = float(bc_value_south.get())

        m_nodes = int(m_input.get())
        n_nodes = int(n_input.get())

        # Generate Non-Dim Numbers
        dx = x_dim/m_nodes
        dy = y_dim/n_nodes

        Aw = dy*z_dim
        Ae = dy*z_dim
        As = dx*z_dim
        An = dx*z_dim
     
        # Generate Coefficients Matrix
        coeff_matrix = np.zeros((7,9))

        # Define Coefficients for Interior Nodes
        coeff_matrix[0,0] = 0.
        coeff_matrix[1,0] = 0.
        coeff_matrix[2,0] = k_const*Aw/dx
        coeff_matrix[3,0] = k_const*Ae/dx
        coeff_matrix[4,0] = k_const*As/dy
        coeff_matrix[5,0] = k_const*An/dy
        coeff_matrix[6,0] = coeff_matrix[2,0]+coeff_matrix[3,0]+coeff_matrix[4,0]+coeff_matrix[5,0]

        # Check West Inputs
        if input_west == 'Conduction':
            coeff_matrix[0,1] = Aw*value_west
            coeff_matrix[1,1] = 0.
            coeff_matrix[2,1] = 0.
            coeff_matrix[3,1] = k_const*Ae/dx
            coeff_matrix[4,1] = k_const*As/dy
            coeff_matrix[5,1] = k_const*An/dy
            coeff_matrix[6,1] = coeff_matrix[2,1]+coeff_matrix[3,1]+coeff_matrix[4,1]+coeff_matrix[5,1]-coeff_matrix[1,1]
            pass
        if input_west == 'Temperature':
            coeff_matrix[0,1] = 2*k_const*Aw/dx*value_west
            coeff_matrix[1,1] = -2*k_const*Aw/dx
            coeff_matrix[2,1] = 0
            coeff_matrix[3,1] = k_const*Ae/dx
            coeff_matrix[4,1] = k_const*As/dy
            coeff_matrix[5,1] = k_const*An/dy
            coeff_matrix[6,1] = coeff_matrix[2,1]+coeff_matrix[3,1]+coeff_matrix[4,1]+coeff_matrix[5,1]-coeff_matrix[1,1]
            pass
        if input_west == 'Insulated':
            coeff_matrix[0,1] = 0.
            coeff_matrix[1,1] = 0.
            coeff_matrix[2,1] = 0.
            coeff_matrix[3,1] = k_const*Ae/dx
            coeff_matrix[4,1] = k_const*As/dy
            coeff_matrix[5,1] = k_const*An/dy
            coeff_matrix[6,1] = coeff_matrix[2,1]+coeff_matrix[3,1]+coeff_matrix[4,1]+coeff_matrix[5,1]-coeff_matrix[1,1]
            pass

        # Check North Inputs
        if input_north == 'Conduction':
            coeff_matrix[0,3] = An*value_north
            coeff_matrix[1,3] = 0.
            coeff_matrix[2,3] = k_const*Aw/dx
            coeff_matrix[3,3] = k_const*Ae/dx
            coeff_matrix[4,3] = k_const*As/dy
            coeff_matrix[5,3] = 0.
            coeff_matrix[6,3] = coeff_matrix[2,3]+coeff_matrix[3,3]+coeff_matrix[4,3]+coeff_matrix[5,3]-coeff_matrix[1,3]
            pass
        if input_north == 'Temperature':
            coeff_matrix[0,3] = 2*k_const*An/dy*value_north
            coeff_matrix[1,3] = -2*k_const*An/dy
            coeff_matrix[2,3] = k_const*Aw/dx
            coeff_matrix[3,3] = k_const*Ae/dx
            coeff_matrix[4,3] = k_const*As/dy
            coeff_matrix[5,3] = 0.
            coeff_matrix[6,3] = coeff_matrix[2,3]+coeff_matrix[3,3]+coeff_matrix[4,3]+coeff_matrix[5,3]-coeff_matrix[1,3]
            pass
        if input_north == 'Insulated':
            coeff_matrix[0,3] = 0.
            coeff_matrix[1,3] = 0.
            coeff_matrix[2,3] = k_const*Aw/dx
            coeff_matrix[3,3] = k_const*Ae/dx
            coeff_matrix[4,3] = k_const*As/dy
            coeff_matrix[5,3] = 0.
            coeff_matrix[6,3] = coeff_matrix[2,3]+coeff_matrix[3,3]+coeff_matrix[4,3]+coeff_matrix[5,3]-coeff_matrix[1,3]
            pass

        # Check East Inputs
        if input_east == 'Conduction':
            coeff_matrix[0,5] = Aw*value_west
            coeff_matrix[1,5] = 0.
            coeff_matrix[2,5] = k_const*Aw/dx
            coeff_matrix[3,5] = 0.
            coeff_matrix[4,5] = k_const*As/dy
            coeff_matrix[5,5] = k_const*An/dy
            coeff_matrix[6,5] = coeff_matrix[2,1]+coeff_matrix[3,1]+coeff_matrix[4,1]+coeff_matrix[5,1]-coeff_matrix[1,1]
            pass
        if input_east == 'Temperature':
            coeff_matrix[0,5] = 2*k_const*Ae/dx*value_east
            coeff_matrix[1,5] = -2*k_const*Ae/dx
            coeff_matrix[2,5] = k_const*Aw/dx
            coeff_matrix[3,5] = 0.
            coeff_matrix[4,5] = k_const*As/dy
            coeff_matrix[5,5] = k_const*An/dy
            coeff_matrix[6,5] = coeff_matrix[2,1]+coeff_matrix[3,1]+coeff_matrix[4,1]+coeff_matrix[5,1]-coeff_matrix[1,1]
            pass
        if input_east == 'Insulated':
            coeff_matrix[0,5] = 0.
            coeff_matrix[1,5] = 0.
            coeff_matrix[2,5] = k_const*Aw/dx
            coeff_matrix[3,5] = 0.
            coeff_matrix[4,5] = k_const*As/dy
            coeff_matrix[5,5] = k_const*An/dy
            coeff_matrix[6,5] = coeff_matrix[2,5]+coeff_matrix[3,5]+coeff_matrix[4,5]+coeff_matrix[5,5]-coeff_matrix[1,5]
            pass

        # Check south Inputs
        if input_south == 'Conduction':
            coeff_matrix[0,7] = As*value_south
            coeff_matrix[1,7] = 0.
            coeff_matrix[2,7] = k_const*Aw/dx
            coeff_matrix[3,7] = k_const*Ae/dx
            coeff_matrix[4,7] = 0.
            coeff_matrix[5,7] = k_const*An/dy
            coeff_matrix[6,7] = coeff_matrix[2,7]+coeff_matrix[3,7]+coeff_matrix[4,7]+coeff_matrix[5,7]-coeff_matrix[1,7]
            pass
        if input_south == 'Temperature':
            coeff_matrix[0,7] = 2*k_const*As/dy*value_south
            coeff_matrix[1,7] = -2*k_const*As/dy
            coeff_matrix[2,7] = k_const*Aw/dx
            coeff_matrix[3,7] = k_const*Ae/dx
            coeff_matrix[4,7] = 0.
            coeff_matrix[5,7] = k_const*An/dy
            coeff_matrix[6,7] = coeff_matrix[2,7]+coeff_matrix[3,7]+coeff_matrix[4,7]+coeff_matrix[5,7]-coeff_matrix[1,7]
            pass
        if input_south == 'Insulated':
            coeff_matrix[0,7] = 0.
            coeff_matrix[1,7] = 0.
            coeff_matrix[2,7] = k_const*Aw/dx
            coeff_matrix[3,7] = k_const*Ae/dx
            coeff_matrix[4,7] = 0.
            coeff_matrix[5,7] = k_const*An/dy
            coeff_matrix[6,7] = coeff_matrix[2,7]+coeff_matrix[3,7]+coeff_matrix[4,7]+coeff_matrix[5,7]-coeff_matrix[1,7]
            pass

        # Generate Corner Nodes
        # Generate West-North Corner
        coeff_matrix[0,2] = coeff_matrix[0,1]+coeff_matrix[0,3]
        coeff_matrix[1,2] = coeff_matrix[1,1]+coeff_matrix[1,3]
        coeff_matrix[2,2] = 0.
        coeff_matrix[3,2] = k_const*Ae/dx
        coeff_matrix[4,2] = k_const*As/dy
        coeff_matrix[5,2] = 0.
        coeff_matrix[6,2] = coeff_matrix[2,2]+coeff_matrix[3,2]+coeff_matrix[4,2]+coeff_matrix[5,2]-coeff_matrix[1,2]

        # Generate North-East Corner
        coeff_matrix[0,4] = coeff_matrix[0,3]+coeff_matrix[0,5]
        coeff_matrix[1,4] = coeff_matrix[1,3]+coeff_matrix[1,5]
        coeff_matrix[2,4] = k_const*Aw/dx
        coeff_matrix[3,4] = 0.
        coeff_matrix[4,4] = k_const*As/dy
        coeff_matrix[5,4] = 0.
        coeff_matrix[6,4] = coeff_matrix[2,4]+coeff_matrix[3,4]+coeff_matrix[4,4]+coeff_matrix[5,4]-coeff_matrix[1,4]

        # Generate East-South Corner
        coeff_matrix[0,6] = coeff_matrix[0,5]+coeff_matrix[0,7]
        coeff_matrix[1,6] = coeff_matrix[1,5]+coeff_matrix[1,7]
        coeff_matrix[2,6] = k_const*Aw/dx
        coeff_matrix[3,6] = 0.
        coeff_matrix[4,6] = 0.
        coeff_matrix[5,6] = k_const*An/dy
        coeff_matrix[6,6] = coeff_matrix[2,6]+coeff_matrix[3,6]+coeff_matrix[4,6]+coeff_matrix[5,6]-coeff_matrix[1,6]

        # Generate South-West Corner
        coeff_matrix[0,8] = coeff_matrix[0,7]+coeff_matrix[0,0]
        coeff_matrix[1,8] = coeff_matrix[1,7]+coeff_matrix[1,0]
        coeff_matrix[2,8] = 0.
        coeff_matrix[3,8] = k_const*Ae/dx
        coeff_matrix[4,8] = 0.
        coeff_matrix[5,8] = k_const*An/dy
        coeff_matrix[6,8] = coeff_matrix[2,8]+coeff_matrix[3,8]+coeff_matrix[4,8]+coeff_matrix[5,8]-coeff_matrix[1,8]


        print(coeff_matrix)

    except ValueError:
        pass

def Reset():
    pass

root = Tk()
root.title("2D Conduction")

mainframe = ttk.Frame(root, padding="1 1 1 1")
mainframe.grid(column=0, row=0, sticky=(W))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

X_input = StringVar()
Y_input = StringVar()
Z_input = StringVar()

k_input = StringVar()

m_input = StringVar()
n_input = StringVar()

bc_input_west = StringVar()
bc_input_north = StringVar()
bc_input_east = StringVar()
bc_input_south = StringVar()

bc_value_west = StringVar()
bc_value_north = StringVar()
bc_value_east = StringVar()
bc_value_south = StringVar()

ttk.Label(mainframe, text="Global Variables").grid(column=1, row=1, sticky=W)
ttk.Label(mainframe, text="X Dimension").grid(column=1, row=2, sticky=W)
ttk.Label(mainframe, text="Y Dimension").grid(column=1, row=3, sticky=W)
ttk.Label(mainframe, text="Z Dimension").grid(column=1, row=4, sticky=W)
ttk.Label(mainframe, text="Coefficient (k)").grid(column=1, row=5, sticky=W)

X_entry = ttk.Entry(mainframe, width=4, textvariable=X_input)
X_entry.grid(column=2, row=2, sticky=(W))

Y_entry = ttk.Entry(mainframe, width=4, textvariable=Y_input)
Y_entry.grid(column=2, row=3, sticky=(W))

Z_entry = ttk.Entry(mainframe, width=4, textvariable=Z_input)
Z_entry.grid(column=2, row=4, sticky=(W))

k_entry = ttk.Entry(mainframe, width=4, textvariable=k_input)
k_entry.grid(column=2, row=5, sticky=(W))

ttk.Label(mainframe, text="(m)").grid(column=3, row=2, sticky=W)
ttk.Label(mainframe, text="(m)").grid(column=3, row=3, sticky=W)
ttk.Label(mainframe, text="(m)").grid(column=3, row=4, sticky=W)
ttk.Label(mainframe, text="(W/m^2/K)").grid(column=3, row=5, sticky=W)

ttk.Label(mainframe, text="Boundary Conditions").grid(column=4, row=1, sticky=W)
ttk.Label(mainframe, text="West Wall").grid(column=4, row=2, sticky=W)
ttk.Label(mainframe, text="North Wall").grid(column=4, row=3, sticky=W)
ttk.Label(mainframe, text="East Wall").grid(column=4, row=4, sticky=W)
ttk.Label(mainframe, text="South Wall").grid(column=4, row=5, sticky=W)

bc_west = ttk.Combobox(mainframe, textvariable=bc_input_west)
bc_west['values'] = ('Conduction', 'Temperature', 'Insulated')
bc_west.set('Conduction')
bc_west.grid(column=5, row=2, sticky=W)

bc_north = ttk.Combobox(mainframe, textvariable=bc_input_north)
bc_north['values'] = ('Conduction', 'Temperature', 'Insulated')
bc_north.set('Conduction')
bc_north.grid(column=5, row=3, sticky=W)

bc_east = ttk.Combobox(mainframe, textvariable=bc_input_east)
bc_east['values'] = ('Conduction', 'Temperature', 'Insulated')
bc_east.set('Conduction')
bc_east.grid(column=5, row=4, sticky=W)

bc_south = ttk.Combobox(mainframe, textvariable=bc_input_south)
bc_south['values'] = ('Conduction', 'Temperature', 'Insulated')
bc_south.set('Conduction')
bc_south.grid(column=5, row=5, sticky=W)

bc_entry_west = ttk.Entry(mainframe, width=4, textvariable=bc_value_west)
bc_entry_west.grid(column=6, row=2, sticky=(W))

bc_entry_north = ttk.Entry(mainframe, width=4, textvariable=bc_value_north)
bc_entry_north.grid(column=6, row=3, sticky=(W))

bc_entry_east = ttk.Entry(mainframe, width=4, textvariable=bc_value_east)
bc_entry_east.grid(column=6, row=4, sticky=(W))

bc_entry_south = ttk.Entry(mainframe, width=4, textvariable=bc_value_south)
bc_entry_south.grid(column=6, row=5, sticky=(W))

ttk.Label(mainframe, text="(W)").grid(column=7, row=2, sticky=W)
ttk.Label(mainframe, text="(C)").grid(column=7, row=3, sticky=W)
ttk.Label(mainframe, text="(na)").grid(column=7, row=4, sticky=W)
ttk.Label(mainframe, text="(na)").grid(column=7, row=5, sticky=W)

ttk.Label(mainframe, text="Solution").grid(column=1, row=15, sticky=W)
ttk.Label(mainframe, text="# of Nodes (N to S)").grid(column=1, row=16, sticky=W)
ttk.Label(mainframe, text="# of Nodes (W to E)").grid(column=1, row=17, sticky=W)

m_entry = ttk.Entry(mainframe, width=4, textvariable=m_input)
m_entry.grid(column=2, row=16, sticky=(W))

n_entry = ttk.Entry(mainframe, width=4, textvariable=n_input)
n_entry.grid(column=2, row=17, sticky=(W))

start_button = ttk.Button(mainframe, text='Start', command=Start)
start_button.grid(column=4, row=16, sticky=(NE))
reset_button = ttk.Button(mainframe, text='Reset', command=Reset)
reset_button.grid(column=5, row=16, sticky=W)

ttk.Label(mainframe, text="Progress:").grid(column=4, row=17, sticky=E)
complete_bar = ttk.Progressbar(mainframe, orient=HORIZONTAL, length=200, mode='indeterminate')
complete_bar.grid(column=5, row=17, sticky=(NE))

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

root.bind('<Return>', Start)
root.mainloop()