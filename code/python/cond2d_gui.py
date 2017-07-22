from tkinter import *
from tkinter import ttk

def	none():
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
bc_west['values'] = ('Conduction', 'Fixed Temperature', 'Insulated')
bc_west.set('Conduction')
bc_west.grid(column=5, row=2, sticky=W)

bc_north = ttk.Combobox(mainframe, textvariable=bc_input_north)
bc_north['values'] = ('Conduction', 'Fixed Temperature', 'Insulated')
bc_north.set('Conduction')
bc_north.grid(column=5, row=3, sticky=W)

bc_east = ttk.Combobox(mainframe, textvariable=bc_input_east)
bc_east['values'] = ('Conduction', 'Fixed Temperature', 'Insulated')
bc_east.set('Conduction')
bc_east.grid(column=5, row=4, sticky=W)

bc_south = ttk.Combobox(mainframe, textvariable=bc_input_south)
bc_south['values'] = ('Conduction', 'Fixed Temperature', 'Insulated')
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

start_button = ttk.Button(mainframe, text='Start', command=none)
start_button.grid(column=4, row=16, sticky=(NE))
reset_button = ttk.Button(mainframe, text='Reset', command=none)
reset_button.grid(column=5, row=16, sticky=W)

ttk.Label(mainframe, text="Progress:").grid(column=4, row=17, sticky=E)
complete_bar = ttk.Progressbar(mainframe, orient=HORIZONTAL, length=200, mode='indeterminate')
complete_bar.grid(column=5, row=17, sticky=(NE))

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

root.mainloop()