import os
import subprocess
import tkinter as tk
from tkinter import messagebox, filedialog
from tkinter import ttk
from tkinter import scrolledtext
import threading
import re

class FGrS_GUI:
    def __init__(self, root):
        self.root = root
        self.root.geometry("475x650")  # Increased height to accommodate output box
        self.root.title("FGrS: Fast Gravity Synthesis")

        self.create_logo()
        self.create_input_fields()
        self.create_buttons()
        self.create_output_box()

        self.exe_path = ""

    def create_logo(self):
        canvas = tk.Canvas(self.root, width=150, height=45, bg='#2c3e50', highlightthickness=0)
        canvas.place(x=150, y=00)
        canvas.create_text(80, 25, text='FGrS', font=('Times New Roman', 26, 'bold'), fill='white')

    def create_input_fields(self):
        # EGM file
        tk.Label(self.root, text="EGM file:").place(x=25, y=50)
        self.egm_file = tk.StringVar(value="D:/data/SH_model/EGM2008.gfc")
        self.egm_file_entry = tk.Entry(self.root, textvariable=self.egm_file, width=45)
        self.egm_file_entry.place(x=120, y=50)
        tk.Button(self.root, text="Browse", command=self.browse_egm_file).place(x=400, y=48)

        # Nmin and Nmax
        tk.Label(self.root, text="Nmin:").place(x=25, y=75)
        self.nmin = tk.Entry(self.root, width=15)
        self.nmin.insert(0, "0")
        self.nmin.place(x=120, y=75)

        tk.Label(self.root, text="Nmax:").place(x=250, y=75)
        self.nmax = tk.Entry(self.root, width=15)
        self.nmax.insert(0, "360")
        self.nmax.place(x=300, y=75)

        # Normal Dropdown
        tk.Label(self.root, text="Normal field:").place(x=25, y=100)
        normal_options = [
            "1 GRS80", "2 WGS84", "3 EGM2008", "4 NIMA96", "5 WGS72",
            "6 GRS67", "7 GRIM5", "8 JGM3", "9 MOON", "10 MARS", "11 VENUS"
        ]
        self.normal = tk.StringVar()
        self.normal.set(normal_options[0])
        self.normal_combobox = ttk.Combobox(self.root, textvariable=self.normal, values=normal_options, width=42)
        self.normal_combobox.place(x=120, y=100)

        # Output Functional Dropdown
        tk.Label(self.root, text="Functional:").place(x=25, y=125)
        output_functional_options = [
            "1 gravity and gravitational potential: W, V",
            "2 anomalous potential: T = W -- U",
            "3 height anomaly: Z = T/gamma_q",
            "4 geoid undulation: N = T0/gamma_0 + ADCB",
            "5 spherical app gravity disturbance: - dT/dr",
            "6 sph app gravity anomaly: 2/r T-dT/dr",
            "7 height: H",
            "8 CS coef defines output: none",
            "9 DOV (xi , eta)",
            "10 gravity anomaly: norm(g) - norm(gamma_q)",
            "11 gravity: gx, gy, gz (LNOF)",
            "12 gravity disturbance: dgx,dgy,dgz (LNOF)",
            "13 gravity disturbance: norm(g)- norm(gamma)",
            "14 dV/dx, dV/dy, dV/dz (LNOF)",
            "15 dV/dphi, dV/dlam, dV/dr",
            "16 T_xx, T_yy, T_zz, T_xy, T_xz, T_yz",
            "17 T_r, T_rr, T_p, T_pp, T_l, T_ll, T_rp, T_rl, T_pl"
        ]
        self.output_functional = tk.StringVar()
        self.output_functional.set(output_functional_options[0])
        self.output_functional_combobox = ttk.Combobox(self.root, textvariable=self.output_functional, values=output_functional_options, width=42)
        self.output_functional_combobox.place(x=120, y=125)

        # Computation Mode
        tk.Label(self.root, text="Computation Mode:").place(x=25, y=150)
        self.computation_mode = tk.StringVar(value="grid")
        tk.Radiobutton(self.root, text="Grid", variable=self.computation_mode, value="grid", command=self.handle_computation_mode_change).place(x=150, y=150)
        tk.Radiobutton(self.root, text="Scatter", variable=self.computation_mode, value="scatter", command=self.handle_computation_mode_change).place(x=220, y=150)
        tk.Radiobutton(self.root, text="Irregular Grid", variable=self.computation_mode, value="irregular-grid", command=self.handle_computation_mode_change).place(x=300, y=150)

        # Coordinates input
        #self.south = self.create_entry("Min Lat (deg):", 7, "30")
        #self.step_south = self.create_entry("Step Lat (arc second):", 8, "60")
        #self.north = self.create_entry("Max Lat (deg):", 9, "40")
        #self.west = self.create_entry("Min Lon (deg):", 10, "30")
        #self.step_west = self.create_entry("Step Lon (arc second):", 11, "60")
        #self.east = self.create_entry("Max Lon (deg):", 12, "40")
        #self.height = self.create_entry("Height (meter):", 13, "0.000")
        
        tk.Label(self.root, text="Min lat. (deg):").place(x=25, y=175)
        self.south = tk.Entry(self.root, width=13)
        self.south.insert(0, "30")
        self.south.place(x=120, y=175)

        tk.Label(self.root, text="Max lat. (deg):").place(x=225, y=175)
        self.north = tk.Entry(self.root, width=13)
        self.north.insert(0, "60")
        self.north.place(x=313, y=175)
  
        tk.Label(self.root, text="Min lon. (deg):").place(x=25, y=200)
        self.west = tk.Entry(self.root, width=13)
        self.west.insert(0, "40")
        self.west.place(x=120, y=200)

        tk.Label(self.root, text="Max lon. (deg):").place(x=225, y=200)
        self.east = tk.Entry(self.root, width=13)
        self.east.insert(0, "65")
        self.east.place(x=313, y=200)

        
        tk.Label(self.root, text="Step lat. (sec):").place(x=25, y=225)
        self.step_south = tk.Entry(self.root, width=13)
        self.step_south.insert(0, "300")
        self.step_south.place(x=120, y=225)

        tk.Label(self.root, text="Step lon. (sec):").place(x=225, y=225)
        self.step_west = tk.Entry(self.root, width=13)
        self.step_west.insert(0, "300")
        self.step_west.place(x=313, y=225)

        tk.Label(self.root, text="Height (m):").place(x=25, y=250)
        self.height = tk.Entry(self.root, width=45)
        self.height.insert(0, "1000.000")
        self.height.place(x=120, y=250)

        
        # Scatter/Irregular File
        tk.Label(self.root, text="Scatter/igrid File:").place(x=25, y=275)
        self.scatter_file = tk.StringVar(value="input.txt")
        self.scatter_file_entry = tk.Entry(self.root, textvariable=self.scatter_file, width=50, state=tk.DISABLED)
        self.scatter_file_entry.place(x=120, y=275)
        tk.Button(self.root, text="Browse", command=self.browse_scatter_file).place(x=400, y=275)

        # Output File
        #self.output_file = self.create_entry("Output File:", 15, "output.txt")
        tk.Label(self.root, text="Output File:").place(x=25, y=300)
        self.output_file = tk.Entry(self.root, width=45)
        self.output_file.insert(0, "output.txt")
        self.output_file.place(x=120, y=300)

        # Thread Number
        #self.thread_number = self.create_entry("Thread Number:", 16, "8")
        tk.Label(self.root, text="Thread Number:").place(x=25, y=320)
        self.thread_number = tk.Entry(self.root, width=45)
        self.thread_number.insert(0, "6")
        self.thread_number.place(x=120, y=320) 

        tk.Label(self.root, text="Results                                       Progress:").place(x=25, y=375)
        progress_frame = tk.Frame(self.root)
        progress_frame.place(x=250, y=375)

        style = ttk.Style()
        
        style.configure("Custom.Horizontal.TProgressbar",
                troughcolor='gray85',
                troughrelief='flat',
                troughborderwidth=1,
                borderwidth=0,
                pbarrelief='flat',
                thickness=1000)
        self.progress_bar = ttk.Progressbar(progress_frame, length=120, mode='determinate',style="Custom.Horizontal.TProgressbar")
        self.progress_bar.pack(side=tk.LEFT, padx=(0, 10))

        self.progress_label = tk.Label(progress_frame, text="0%")
        self.progress_label.pack(side=tk.LEFT)
 
        self.coord_widgets = [self.south, self.step_south, self.north, self.west, self.step_west, self.east, self.height]

    def create_entry(self, label, row, default_value):
        tk.Label(self.root, text=label).grid(row=row, column=0, sticky='e', pady=2)
        entry = tk.Entry(self.root, width=50)
        entry.insert(0, default_value)
        entry.grid(row=row, column=1)
        return entry

    def create_buttons(self):
        tk.Button(self.root, text="RUN", command=self.run_fgrs, bg="green", fg="white", width=12).place(x=25,y=350)
        tk.Button(self.root, text="EXIT", command=self.root.quit, bg="red", fg="white", width=12).place(x=300,y=350)
        
    def create_output_box(self):
        self.output_box = scrolledtext.ScrolledText(self.root, wrap=tk.WORD, width=53, height=13)
        self.output_box.place(x=25, y=400)

    def browse_egm_file(self):
        egm_file_path = filedialog.askopenfilename(title="Select EGM File", filetypes=[("EGM files", "*.gfc"), ("All files", "*.*")])
        if egm_file_path:
            self.egm_file.set(egm_file_path)

    def browse_scatter_file(self):
        scatter_file_path = filedialog.askopenfilename(title="Select Scatter/Irregular Grid File", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if scatter_file_path:
            self.scatter_file.set(scatter_file_path)

    def handle_computation_mode_change(self):
        if self.computation_mode.get() == "grid":
            self.scatter_file_entry.config(state=tk.DISABLED)
            self.coord_entries_state(tk.NORMAL)
        else:
            self.scatter_file_entry.config(state=tk.NORMAL)
            self.coord_entries_state(tk.DISABLED)

    def coord_entries_state(self, state):
        for widget in self.coord_widgets:
            widget.config(state=state)
 
    
    def run_fgrs(self):
        self.output_box.delete(1.0, tk.END)
        self.progress_bar['value'] = 0
        self.progress_label.config(text="0%")
        
        if not self.validate_inputs():
            return
        
        self.write_ioptions_file()
        self.exe_path = os.path.join(os.getcwd(), "FGrS.exe")
        
        if not self.exe_path:
            self.exe_path = filedialog.askopenfilename(title="Select FGrS.exe", filetypes=(("Executable Files", "*.exe"),))

        if not self.exe_path:
            self.output_box.insert(tk.END, "No executable file selected.\n")
            return
        
        try:
            process = subprocess.Popen([self.exe_path, "ioptions.txt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, creationflags=subprocess.CREATE_NO_WINDOW)

            def read_output():
                while True:
                    output = process.stdout.readline()
                    clean_output = re.sub(r'\x08', '', output) 
                    if clean_output.startswith("progress"):
                        if clean_output.startswith("progress completed"):
                            progress = 100
                            self.update_progress(progress)
                            clean_output="progress completed"
                            self.output_box.insert(tk.END, clean_output)
                            self.output_box.see(tk.END)
                            continue
                        else:
                            progress = int(clean_output.split()[1])
                        self.update_progress(progress)
                        continue
                    if clean_output == '' and process.poll() is not None:
                        break
                    if clean_output:
                        self.output_box.insert(tk.END, clean_output)
                        self.output_box.see(tk.END)
                process.stdout.close()

            threading.Thread(target=read_output, daemon=True).start()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")
            self.output_box.insert(tk.END, f"Error: {e}\n")
            self.output_box.see(tk.END)

    def update_progress(self, progress):
        self.progress_bar['value'] = progress
        self.progress_label.config(text=f"{progress}%")
        self.root.update_idletasks()

    def validate_inputs(self):
        try:
            nmin_val = int(self.nmin.get())
            nmax_val = int(self.nmax.get())
            if nmin_val < 0 or nmax_val < 0 or nmin_val >= nmax_val:
                raise ValueError("nmin must be a non-negative integer less than nmax.")

            south_val = float(self.south.get())
            north_val = float(self.north.get())
            west_val = float(self.west.get())
            east_val = float(self.east.get())

            if not (-90 < south_val < 90) or not (-90 < north_val < 90):
                raise ValueError("Latitude must be in the range [-90, 90].")
            if not (-180 <= west_val <= 360) or not (-180 <= east_val <= 360):
                raise ValueError("Longitude must be in the range [-180, 180] or [0, 360].")
            if south_val >= north_val:
                raise ValueError("Min latitude must be less than max latitude.")
            if west_val >= east_val:
                raise ValueError("Min longitude must be less than max longitude.")

            float(self.height.get())
            int(self.thread_number.get())

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return False
        return True

    def write_ioptions_file(self):
        with open("ioptions.txt", "w") as f:
            f.write(f"{self.egm_file.get()}\n")
            f.write(f"{self.nmin.get()}\n")
            f.write(f"{self.nmax.get()}\n")
            f.write(f"{self.normal.get().split()[0]}\n")
            f.write(f"{self.output_functional.get().split()[0]}\n")
            f.write(f"{self.computation_mode.get()}\n")
            f.write(f"{self.south.get()}\n")
            f.write(f"{self.step_south.get()}\n")
            f.write(f"{self.north.get()}\n")
            f.write(f"{self.west.get()}\n")
            f.write(f"{self.step_west.get()}\n")
            f.write(f"{self.east.get()}\n")
            f.write(f"{self.height.get()}\n")
            f.write(f"{self.scatter_file.get()}\n")
            f.write(f"{self.output_file.get()}\n")
            f.write(f"{self.thread_number.get()}\n")
            f.close()

if __name__ == "__main__":
    root = tk.Tk()
    gui = FGrS_GUI(root)
    root.mainloop()
