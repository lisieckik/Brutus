import time
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from requests import request
import numpy as np
import h5py, gc, os
from astropy.io import fits
from astropy.table import Table

# Function for opening the
# file explorer window
def browseFiles(event):
    filename = filedialog.askopenfilename(initialdir=os.getcwd(),
                                          title="Select a File",)
    if len(filename) != 0:
        if event.widget == button5_browse:
            text_new = open(filename, 'r').read()
            para_box.delete("1.0","end")
            para_box.insert('end', text_new)
            check_text()
        elif event.widget == load:
            global id_file_load
            id_file_load = filename
            recover(filename)
        try:
            if event.widget == file_id:
                ids_text.delete(0, tk.END)
                ids_text.insert('end', filename)
                global_id_file(file_id)
        except:pass
    return 'break'
def recover(oldone):
    global id_file_load
    try:
        oldone = open(oldone, 'r').read().split('###\n')
        try:
            input_cats_box.delete(0, tk.END)
        except:
            pass
        input_cats_box.insert('end', oldone[0][0:-1])

        try:
            para_box.delete("1.0", "end")
        except:
            pass
        paras_load = oldone[1].split('\n')
        for i in paras_load:
            if len(i) > 0:
                para_box.insert('end', i + '\n')
        check_text()

        try:
            e1.delete(0, tk.END)
        except:
            pass
        try:
            e2.delete(0, tk.END)
        except:
            pass
        try:
            e3.delete(0, tk.END)
        except:
            pass
        try:
            e4.delete(0, tk.END)
        except:
            pass

        e1.insert('end', oldone[2].split()[0])
        e2.insert('end', oldone[2].split()[1])
        e3.insert('end', oldone[2].split()[2])
        e4.insert('end', oldone[2].split()[3])

        try:
            out_cats_text.delete(0, tk.END)
        except:
            pass
        out_cats_text.insert('end', oldone[3][0:-1])

        ids_load = oldone[5].split()
        try:
            if len(ids_load[1]) > 0:
                id_true.set(int(ids_load[0]))
                if fake2.winfo_exists() == 1:
                    id_file_load = ids_load[1]
                    select_id()
        except:
            pass
    except:pass
def intput_cat(event):
    filename = filedialog.askdirectory(initialdir=os.getcwd(),
                                        title="Select a Directory",)
    if len(filename) != 0:
        filename += '/'
        if event.widget == button1_explore:
            print(1)
            input_cats_box.delete(0, tk.END)
            input_cats_box.insert(0,filename)
        elif event.widget == button_out:
            out_cats_text.delete(0, tk.END)
            out_cats_text.insert(0,filename)
    return 'break'
def save_new_file(event):
    if event.widget == save:
        f = filedialog.asksaveasfile(initialdir=os.getcwd(),initialfile='Untitled.txt',
                      defaultextension=".txt", filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt")])
        if f != None:
            write_before_saving(f)
    return 'break'
def save_when_close():
    f = open('last_used_config.txt', 'w')
    write_before_saving(f)
    root.destroy()
def write_before_saving(f):
    f.write(input_cats_box.get())
    f.write('\n###\n')
    f.write(para_box.get('1.0', 'end')[0:-1])
    f.write('\n###\n')
    f.write('%s %s %s %s' % (e1.get(), e2.get(), e3.get(), e4.get()))
    f.write('\n###\n')
    f.write(out_cats_text.get())
    f.write('\n###\n')
    if prog_var.get() == 0:
        f.write('0')
    elif prog_var.get() == 1:
        f.write('0')
    f.write('\n###\n')
    if id_true.get() == 0:
        f.write('0\n')
        f.write(ids_text.get())
    elif id_true.get() == 1:
        f.write('1')
def prepare_one_row(i, name, root):
    ti = tk.Label(root, text=name)
    ti.grid(row=i, column=1, sticky='nw', padx=10, pady=1, columnspan=10)
    cbi = tk.Checkbutton(root, width=1, variable=checkbox_value[i], onvalue=1, offvalue=0)
    cbi.grid(row=i, column=0, sticky='nw', padx=10, pady=1, columnspan=1)
    checkbox_list[i]=cbi
def openNewWindow():
    global checkbox_value, names, checkbox_list, newWindow, browser_open
    if browser_open:
        pass
    else:
        browser_open = True
        newWindow = tk.Toplevel(root)
        newWindow.protocol('WM_DELETE_WINDOW', end_chosing)

        # sets the title of the
        # Toplevel widget
        newWindow.title("Browse SIMBA quantities")

        # sets the geometry of toplevel
        newWindow.geometry("500x700")

        names = open("parameters_all.txt").readlines()
        names = [' '.join(k.split()) for k in names]

        ## create a main frame
        left = tk.Frame(newWindow)
        left.pack(fill=tk.BOTH, expand = 1)

        ## create a canvas
        my_canvas = tk.Canvas(left)
        my_canvas.pack(side=tk.LEFT, fill = tk.BOTH, expand=1)

        ## add a scrollbar
        scroll_v = ttk.Scrollbar(left, orient=tk.VERTICAL, command=my_canvas.yview)
        scroll_v.pack(side=tk.RIGHT, fill=tk.Y)

        ## configure the canvas
        my_canvas.configure(yscrollcommand=scroll_v.set)
        my_canvas.bind('<Configure>', lambda e: my_canvas.configure(
            scrollregion=my_canvas.bbox('all')))

        ## create another frame inside the canvas
        right = tk.Frame(my_canvas)

        ## add that to the canvas
        my_canvas.create_window((0,0), window=right, anchor='nw')

        try:
            checkbox_value.keys()
            for i in range(len(names)):
                try: checkbox_value[i]
                except:
                    checkbox_value[i] = tk.IntVar(root, 0)
        except:
            checkbox_value = {}
            checkbox_list = {}
            for i in range(len(names)):
                checkbox_value[i] = tk.IntVar(root, 0)
        for i in range(len(names)):
            prepare_one_row(i, names[i], right)

        save = tk.Button(newWindow, text='Save', command=end_chosing)
        save.pack(anchor='se',padx=30, pady=10, fill=tk.X)

def end_chosing():
    global browser_open
    text = ''
    for i in range(len(checkbox_value)):
        if checkbox_value[i].get():
            text += names[i] + '\n'
    para_box.delete("1.0","end")
    para_box.insert('end', text)
    browser_open = False
    newWindow.destroy()
def check_text():
    global checkbox_value, checkbox_list
    try: checkbox_value.keys()
    except:
        checkbox_value = {}
        checkbox_list = {}
    tt = para_box.get('1.0', 'end').split('\n')[0:-1]
    names = open("parameters_all.txt").readlines()
    names = [' '.join(k.split()) for k in names]
    wrong = []
    for i in range(len(tt)):
        found = False
        words = tt[i].split()
        if len(words)==0:
            found = True
        else:
            for j in range(len(names)):
                words_n = names[j].split()
                if len(words) == 2 and len(words_n) ==2:
                    if words_n[0] == words[0] and words_n[1] == words[1]:
                        found = True
                        checkbox_value[j] = tk.IntVar(root, 1)
                        break
                elif len(words) == 3 and len(words_n) ==3:
                    if words_n[0] == words[0] and words_n[1] == words[1] and words_n[2] == words[2]:
                        found = True
                        checkbox_value[j] = tk.IntVar(root, 1)
                        break
        if not found:
            print(i, tt[i])
            wrong.append(i)

    text = para_box.get('1.0', 'end').split('\n')[0:-1]
    text_new = ''
    for k in range(len(text)):
        good = True
        for i in wrong:
            if k == i:
                good = False
                break
        if good:
            if len(text[k].split()) >0:
                text_new += ' '.join(text[k].split()) + '\n'
    if len(text_new) !=0:
        while text_new[-1] == '\n': text_new = text_new[0:-1]
    para_box.delete("1.0","end")
    para_box.insert('end', text_new)

    if len(wrong) !=0:
        newWindow = tk.Toplevel(root)
        newWindow.title("Warning!")

        # sets the geometry of toplevel
        newWindow.geometry("250x250")
        war = tk.Label(newWindow, text ='Following attributes\nare not in SIMBA:')
        war.grid(row = 0, column = 0,sticky='news',padx = 10,
                     pady = 10, columnspan=2)

        scroll_3 = tk.Scrollbar(newWindow, orient='vertical', width=30)
        bad_values = tk.Text(newWindow, width=25, yscrollcommand=scroll_3.set, height=10)
        bad_values.grid(row=2, column=0, sticky='ew', padx=10,
                      pady=0, columnspan=1, rowspan=1)
        scroll_3.grid(row=2, column=1, sticky='ns', padx=0)
        scroll_3.config(command=bad_values.yview)
        text = ''
        for i in wrong:
            text += tt[i] + '\n'
        text = text[0:-1]
        bad_values.insert('end', text)
def onClick(event):
    global entry_number
    if event.widget == e1:
        entry_number = 1
    elif event.widget == e2:
        entry_number = 2
    elif event.widget == e3:
        entry_number = 3
    elif event.widget == e4:
        entry_number = 4
def update_entry():
    global warning_prog
    try:
        if entry_number == 1:
            e3.delete(0,tk.END)
            e3.insert(0,round(z_snap[int(e1.get())],3))
        elif entry_number == 2:
            e4.delete(0,tk.END)
            e4.insert(0,round(z_snap[int(e2.get())],3))
        elif entry_number == 3:
            z_dist = np.abs(z_snap - float(e3.get()))
            z_dist = int(np.argwhere(z_dist==np.min(z_dist))[0])
            e1.delete(0,tk.END)
            e1.insert(0,z_dist)
        elif entry_number == 4:
            z_dist = np.abs(z_snap - float(e4.get()))
            z_dist = int(np.argwhere(z_dist==np.min(z_dist))[0])
            e2.delete(0,tk.END)
            e2.insert(0,z_dist)
    except: pass
    try:
        if prog_var.get() == 0 and id_true.get() == 0 and warning_prog == False:
            warning_prog = True
            newWindow = tk.Toplevel(root)
            newWindow.title("Warning!")
            newWindow.geometry("300x150")
            war = tk.Label(newWindow, text='Be careful!\nIf you specify IDs\n and run it without tracking\nyou can encounter a problem.\n'
                                           'If you have lists of IDs\nfor more than one snapshot\nyou have to run Brutus\none by one.')
            war.pack(fill=tk.BOTH, expand=1)
    except:pass
    entry_frame.after(1000, update_entry)
def global_id_file(event):
    global id_file_load
    try:
        id_file_load = ids_text.get()
    except: pass
def select_id():
    global ids_text, scroll_h4, label_id,file_id, fake1, fake2,prog_var
    if id_true.get() == 0:
        label_id = tk.Label(frame1, text='Select file containing IDs (.txt with one column only)')
        label_id.grid(row=14, column=0, sticky='ew', padx=0, pady=2, columnspan=7)
        scroll_h4 = tk.Scrollbar(frame1, orient='horizontal', width=25)
        ids_text = tk.Entry(frame1, width=50, xscrollcommand=scroll_h4.set)
        ids_text.grid(row=15, column=0, sticky='ew', padx=20, pady=2, columnspan=7)
        ids_text.bind('<FocusOut>', global_id_file)
        scroll_h4.grid(row=16, column=0, columnspan=7, sticky='news', padx=15,pady=0)
        scroll_h4.config(command=ids_text.xview)
        file_id = tk.Button(frame1, text='File explorer')
        file_id.grid(row=14, column= 9, sticky='new', padx = 10, pady = 0, columnspan=2)
        file_id.bind('<Button-1>', browseFiles)
        fake1.destroy()
        fake2.destroy()
        prog_var.set(1)
        try:
            ids_text.insert(0,id_file_load)
        except:pass
    else:
        ids_text.destroy()
        scroll_h4.destroy()
        label_id.destroy()
        file_id.destroy()
        fake1 = tk.Label(frame1, text='\n')
        fake1.grid(row=14, column=9, sticky='news', padx=12,
                   pady=0, columnspan=2, )
        fake2 = tk.Label(frame1, text='\n')
        fake2.grid(row=15, column=9, sticky='news', padx=12,
                   pady=1, columnspan=2, )
def writting_off(event):
    check_text()
def call_brutus():
    warnings = ''
    try:
        files = os.listdir(input_cats_box.get())
        no_hdf5 = True
        for i in range(len(files)):
            if len(files[i].split('.'))>1:
                if files[i].split('.')[-1]=='hdf5':
                    no_hdf5=False
        if no_hdf5:
            warnings += 'There are no hdf5 files in input directory!\n'
    except:
        warnings += 'No input directory specified!\n'

    if len(para_box.get("1.0",tk.END)) ==1:
        warnings += 'No quantities specified!\n'
    else:check_text()

    try:
        int(e1.get())
    except:
        warnings += 'No starting snapshot specified!\n'

    try:
        int(e2.get())
    except:
        warnings += 'No ending snapshot specified!\n'

    try:
        out_cats_text.get()
    except:
        warnings += 'Wrong output catalog!\n'

    try:
        input_ids = ids_text.get()
        if input_ids.split('.')[1] == 'txt':
            pass
        else:
            warnings += 'Wrong IDs input extension!\n'
    except:
        pass

    if len(warnings)!=0:
        newWindow = tk.Toplevel(root)
        newWindow.title("Warning!")

        # sets the geometry of toplevel
        newWindow.geometry("375x150")
        war = tk.Label(newWindow, text ='Following problems with BRUTUS input were found:')
        war.grid(row = 0, column = 0,sticky='news',padx = 10,
                     pady = 10, columnspan=2)

        scroll_8 = tk.Scrollbar(newWindow, orient='vertical', width=30)
        bad_values = tk.Text(newWindow, width=40, yscrollcommand=scroll_8.set, height=6)
        bad_values.grid(row=2, column=0, sticky='ew', padx=10,
                      pady=0, columnspan=1, rowspan=1)
        scroll_8.grid(row=2, column=1, sticky='ns', padx=0)
        scroll_8.config(command=bad_values.yview)
        bad_values.insert('end', warnings)
    else:
        try:
            input_ids = ids_text.get()
        except:
            input_ids=-1
        real_brutus(input_cats_box.get(),
                    para_box.get("1.0",tk.END),
                    input_ids,
                    prog_var.get(),
                    [int(e1.get()), int(e2.get())],
                    out_cats_text.get())
    return 'break'
def real_brutus(SIMBA_catalogs_path,
                used_parameters,
                interesting_ids,
                progenitor_track,
                snapshot_range,
                output):

    if progenitor_track == 0: progenitor_track = False
    elif progenitor_track ==1: progenitor_track = True
    used_parameters = used_parameters.split('\n')
    while used_parameters[-1] == '': used_parameters = used_parameters[0:-1]

    ### load IDs of galaxies at the lowest redshift on interest
    if interesting_ids != -1:
        interesting_ids = open(interesting_ids, 'r').read().split()
        slices = np.array(interesting_ids, dtype=int)
        ind_valid = np.where(slices >= 0)[0]
        slices = slices[ind_valid]
        del ind_valid
    ### make output folder if not present
    output_f = os.listdir(output)
    new_f = True
    for f in output_f:
        if f == 'Brutus_output':
            new_f = False
            break
    if new_f: os.mkdir(output + '/Brutus_output')
    del output_f, new_f

    ### retrive attributes for new catalogs
    names = ['redshift']
    attributes_galaxies = []
    halo_index = False


    columns = used_parameters
    for c in columns:
        words = c.split()
        if words[0] == 'halo_data' and not halo_index:
            for c1 in columns:
                w1 = c1.split()
                if w1[1] == 'parent_halo_index': halo_index = True
                break
            if not halo_index:
                names.append('galaxy_data_parent_halo_index')
                attributes_galaxies.append(['galaxy_data', 'parent_halo_index'])
                halo_index = True
            del c1, w1
        name = ''
        attribute = []
        for w in words:
            attribute.append(w)
            name += w + '_'
            if w == 'pos':
                name = 'pos'
            elif w == 'vel':
                name = 'vel'
        attributes_galaxies.append(attribute)
        if name == 'pos':
            names.append(words[0] + '_' + 'pos_x')
            names.append(words[0] + '_' + 'pos_y')
            names.append(words[0] + '_' + 'pos_z')
        elif name == 'vel':
            names.append(words[0] + '_' + 'vel_x')
            names.append(words[0] + '_' + 'vel_y')
            names.append(words[0] + '_' + 'vel_z')
        else:
            name = name[0:-1]
            names.append(name)
    del name, attribute, c, w, columns, halo_index
    if progenitor_track:
        if interesting_ids != -1:
            if redshift_range[0]:
                names.append('descendant_galaxy_at_z_%.3f' % redshift_range[1])
            elif not redshift_range[0]:
                names.append('descendant_galaxy_at_z_%.3f' % z_snap[snapshot_range[1]])
        else:
            names.append('progenitor_galaxy_at_previous_snapshot')
    number_of_columns = len(names)

    i = 0
    j = 0
    while True:
        if i >= 152: break

        ### preparing for chosing the catalogs

        ind_snap = 151 - i
        redshift = z_snap[ind_snap]
        ### checking if the catalog i should be used
        work = False
        if snapshot_range[0] <= snapshot_range[1]:
            if ind_snap >= snapshot_range[0] and ind_snap <= snapshot_range[1]:
                work = True
        else:
            if ind_snap <= snapshot_range[0] and ind_snap >= snapshot_range[1]:
                work = True
        ### if the calatog is used
        if work:
            ind_snap = str(ind_snap)
            if len(ind_snap) < 3:
                ind_snap = '0' * (3 - len(
                    ind_snap)) + ind_snap  # if snapshot number is less than 100 we need to add 0 in front
            status.delete(0, tk.END)
            status.insert('end', 'Now proccesing %s'%ind_snap)
            ### open catalog
            file = SIMBA_catalogs_path + "m100n1024_%s.hdf5" % ind_snap
            with h5py.File(file) as sim:

                ### preparing tables to handle IDs and progenitor-ID relation
                if interesting_ids == -1:
                    slices = np.arange(0, len(sim['galaxy_data']['GroupID']), dtype=int)
                    if progenitor_track:
                        prog_desc = np.empty([len(sim['galaxy_data']['GroupID']), 2], dtype=int)
                        prog_desc[:, 0] = np.arange(0, len(sim['galaxy_data']['GroupID']), dtype=int)
                        prog_desc[:, 1] = np.arange(0, len(sim['galaxy_data']['GroupID']), dtype=int)
                else:
                    if progenitor_track and j == 0:
                        prog_desc = np.empty([len(slices), 2], dtype=int)
                        prog_desc[:, 0] = slices
                        prog_desc[:, 1] = slices

                ### prepare the output table
                if interesting_ids != -1:
                    hdu_out = np.zeros([len(slices), number_of_columns])
                else:
                    hdu_out = np.zeros([len(sim['galaxy_data']['GroupID']), number_of_columns])

                if progenitor_track:
                    if interesting_ids == -1:
                        progenitor_list = sim['tree_data']['progen_galaxy_star'][:, 0]
                        prog_desc[:, 0] = progenitor_list
                    else:
                        progenitor_list = np.array(sim['tree_data']['progen_galaxy_star'])[slices][:, 0]
                        prog_desc[:, 0] = progenitor_list

                ### input redshift
                hdu_out[:, 0] = [z_snap[int(ind_snap)]] * hdu_out.shape[0]

                ### input other quantities
                n = 1
                for k in attributes_galaxies:
                    ### as the slices has to sorted, both, before and after sorting positions, are stored
                    # to keep halo ids with the correct galaxy
                    if k[0] == 'halo_data':
                        subsample = parent_halo_index
                        sorted_subsample = np.argsort(subsample)
                        subsample = subsample[sorted_subsample]
                    else:
                        if interesting_ids == -1:
                            subsample = np.arange(0, len(sim['galaxy_data']['GroupID']), dtype=int)
                        else:
                            subsample = slices
                        sorted_subsample = np.argsort(subsample)
                        subsample = subsample[sorted_subsample]

                    ### inputing quantities according to the attributes
                    if len(k) == 2:
                        sub_data = np.array(sim[k[0]][k[1]])[subsample]
                        if k[1] == 'parent_halo_index':
                            parent_halo_index = sub_data

                        if k[1] == 'pos' or k[1] == 'vel':
                            hdu_out[sorted_subsample, n:n + 3] = sub_data
                            n += 3
                        else:
                            hdu_out[sorted_subsample, n] = sub_data
                            n += 1
                    if len(k) == 3:
                        if interesting_ids == -1:
                            sub_data = sim[k[0]][k[1]][k[2]]
                        else:
                            sub_data = np.array(sim[k[0]][k[1]][k[2]])[subsample]
                        hdu_out[sorted_subsample, n] = sub_data
                        n += 1

                if progenitor_track:
                    hdu_out[:, n] = prog_desc[:, 1]
                    ind_valid = np.where(np.array(prog_desc[:, 0]) != -1)[0]
                    prog_desc = prog_desc[ind_valid]
                    slices = np.array(prog_desc[:, 0])

            cols = Table(hdu_out, names=names)
            cols.write(output + '/Brutus_output/%s_resized.fits' % ind_snap, format='fits', overwrite=True)
            del cols
            j += 1
        i += 1
    status.delete(0, tk.END)
    status.insert('end', 'Finished!')
    status.configure(bg='#0a6e24')

# Initiate the GUI
root = tk.Tk()
root.title('Brutus')
root.geometry("800x850")
root.protocol('WM_DELETE_WINDOW', save_when_close)
# Left frame (basic entry and rows of widgets)

# ## create a main frame
main_frame = tk.Frame(root)
main_frame.pack(fill=tk.BOTH, expand = 1)
## create a canvas
main_canvas = tk.Canvas(main_frame)
main_canvas.pack(side=tk.LEFT, fill = tk.BOTH, expand=1)
## add a scrollbar
scroll_main_v = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=main_canvas.yview)
scroll_main_v.pack(side=tk.RIGHT, fill=tk.Y)

## configure the canvas
main_canvas.configure(yscrollcommand=scroll_main_v.set)
main_canvas.bind('<Configure>', lambda e: main_canvas.configure(
    scrollregion=main_canvas.bbox('all')))
## create another frame inside the canvas
frame1 = tk.Frame(main_canvas)
## add that to the canvas
main_canvas.create_window((0,0), window=frame1, anchor='nw')

### inpit catalogs
input_cats_text = tk.Label(frame1, text='Where are original SIMBA catalogues stored?')
input_cats_text.grid(row = 0, column = 0, sticky='news',padx = 10,
                     pady = 0, columnspan=7)

scroll_h=tk.Scrollbar(frame1, orient='horizontal', width=20)
input_cats_box = tk.Entry(frame1, width=50, xscrollcommand=scroll_h.set)
input_cats_box.grid(row=1, column=0,  sticky='ewn', padx = 20, pady = 10, columnspan=7)
scroll_h.grid(row=2, column=0,columnspan=7,sticky='new', padx=15)
scroll_h.config(command=input_cats_box.xview)

tutorial0 = tk.Label(frame1, text='Download here: http://simba.roe.ac.uk/')
tutorial0.grid(row = 1, column = 9, sticky='ew',padx = 10,
                     pady = 0, rowspan=1, columnspan=2)

browser_open = False
button1_explore = tk.Button(frame1, text='File explorer')
button1_explore.grid(row=0, column= 9, sticky='news', padx = 10, pady = 10, columnspan=2)
button1_explore.bind('<Button-1>', intput_cat)


### used parameters
para_text = tk.Label(frame1, text='What are the quantities of your interest?')
para_text.grid(row = 3, column = 0, sticky='news',padx = 10,
                     pady = 10, columnspan=7)
button2_browse = tk.Button(frame1, text='Browse quantities', command= openNewWindow)
button2_browse.grid(row=3, column= 9, sticky='news', padx = 10, pady = 10, columnspan=1)

button5_browse = tk.Button(frame1, text='From file', command= openNewWindow)
button5_browse.grid(row=3, column= 10, sticky='news', padx = 10, pady = 10, columnspan=1)
button5_browse.bind('<Button-1>', browseFiles)

scroll_2=tk.Scrollbar(frame1, orient='vertical',width=30)
para_box = tk.Text(frame1, width=50, yscrollcommand=scroll_2.set, height= 10)
para_box.grid(row=4, column=0,  sticky='ew', padx = 20,
              pady = 0, columnspan=7, rowspan=1)
para_box.bind('<FocusOut>', writting_off)

scroll_2.grid(row=4, column=8,sticky='ns',padx=0)
scroll_2.config(command=para_box.yview)

tutorial = tk.Label(frame1, text='You can copy and paste\nrows from file '
                                 '\"parameters_all.txt\"\nor select from the browser\n\n'
                                 'Rows have 2 or 3 elemets\n'
                                 'They have to be separated by space\n\n'
                                 'For description check\n'
                                 'CAEZAR documentation')
tutorial.grid(row = 4, column = 9, sticky='ew',padx = 10,
                     pady = 10, rowspan=1, columnspan=2)

button3_browse = tk.Button(frame1, text='Check and save', command= check_text)
button3_browse.grid(row=5, column= 0, sticky='news', padx = 20, pady = 10, columnspan=8)


## snapshot range

### load redshift information
response = request("GET", 'http://simba.roe.ac.uk/outputs.txt')
data = response.text.split('\n')[0:-1]
data = [list(map(float, i.split())) for i in data]
z_snap = np.array(data)[:,1]
response.close()
del response, data

range_text = tk.Label(frame1, text='Which snapshots do you want to study?')
range_text.grid(row = 6, column = 0, sticky='n',padx = 10,
                     pady = 0, columnspan=7)

entry_frame = tk.Frame(frame1)

range_text1 = tk.Label(entry_frame, text='From')
range_text1.grid(row = 0, column = 1, sticky='nw',padx = 60,
                     pady = 4, columnspan=1)

range_text2 = tk.Label(entry_frame, text='To')
range_text2.grid(row = 0, column = 2, sticky='nw',padx = 60,
                     pady = 4, columnspan=1)


range_text3 = tk.Label(entry_frame, text='Snapshot')
range_text3.grid(row = 1, column = 0, sticky='nw',padx = 10,
                     pady = 4, columnspan=1)

range_text4 = tk.Label(entry_frame, text='Redshift')
range_text4.grid(row = 2, column = 0, sticky='nw',padx = 10,
                     pady = 4, columnspan=1)
e1 = tk.Entry(entry_frame, width=15)
e1.grid(row = 1, column = 1,pady=4,padx=7, sticky='n')
e1.bind('<Button-1>', onClick)

e2 = tk.Entry(entry_frame, width=15)
e2.grid(row = 1, column = 2,pady=4,padx=7, sticky='n')
e2.bind('<Button-1>', onClick)

e3 = tk.Entry(entry_frame, width=15)
e3.grid(row = 2, column = 1,pady=4,padx=7, sticky='n')
e3.bind('<Button-1>', onClick)

e4 = tk.Entry(entry_frame, width=15)
e4.grid(row = 2, column = 2,pady=4,padx=7, sticky='n')
e4.bind('<Button-1>', onClick)

entry_frame.grid(row = 7, column = 0, sticky='news',padx = 10,
                     pady = 4, columnspan=7, rowspan=3)
update_entry()
tutorial2 = tk.Label(frame1, text='If you insert snapshots\nthe exact redshift will be found\n\n'
                                  'If you insert redshift\nthe closest snapshot will be found')
tutorial2.grid(row = 7, column = 9, sticky='ew',padx = 10,
                     pady = 0, rowspan=3, columnspan=2)

## output
out_cats_text = tk.Label(frame1, text='Where to store result? New directory will be build.')
out_cats_text.grid(row = 10, column = 0, sticky='news',padx = 10,
                     pady = 10, columnspan=7)
button_out = tk.Button(frame1, text='File explorer')
button_out.grid(row=10, column= 9, sticky='news', padx = 10, pady = 10, columnspan=2)
button_out.bind('<Button-1>', intput_cat)

scroll_h3=tk.Scrollbar(frame1, orient='horizontal', width=25)
out_cats_text = tk.Entry(frame1, width=50, xscrollcommand=scroll_h3.set)
out_cats_text.grid(row=11, column=0,  sticky='ew', padx = 20, pady = 0, columnspan=7)
out_cats_text.insert(0, os.getcwd())
scroll_h3.grid(row=12, column=0,columnspan=7,sticky='news',padx=15)
scroll_h3.config(command=out_cats_text.xview)



## Progenitor track
warning_prog = False
prog = tk.Frame(frame1)
prog_text = tk.Label(prog, text='Do you want to track progenitors?')
prog_text.grid(row = 0, column = 1, sticky='ne',padx = 10,
                     pady = 10, columnspan=2)
prog_var = tk.IntVar()
prog_var.set(0)
prog_box = tk.Checkbutton(prog, width=1, onvalue=1, offvalue=0, variable=prog_var)
prog_box.grid(row = 0, column = 0, sticky='nw',padx = 10,
                     pady = 10, columnspan=1)

id_text = tk.Label(prog, text='Do you want to Study all galaxies?')
id_text.grid(row = 1, column = 1, sticky='ne',padx = 10,
                     pady = 10, columnspan=2)
id_true = tk.IntVar()
id_true.set(1)

id_var = tk.Checkbutton(prog, width=1, variable=id_true, onvalue=1, offvalue=0, command=select_id)
id_var.grid(row = 1, column = 0, sticky='nw',padx = 10,
                     pady = 10, columnspan=1)

prog.grid(row = 13, column = 0, sticky='nw',padx = 40,
                     pady = 0, columnspan=7,rowspan=1)


## finish
fake1 = tk.Label(frame1, text='\n')
fake1.grid(row = 14, column = 9, sticky='news',padx = 12,
                     pady = 1, columnspan=2,)
fake2 = tk.Label(frame1, text='\n')
fake2.grid(row = 15, column = 9, sticky='news',padx = 12,
                     pady = 0, columnspan=2,)


save = tk.Button(frame1, text='Save config')
save.grid(row = 17, column = 9, sticky='nesw',padx = 10,
                     pady = 0, columnspan=1)
save.bind('<Button-1>', save_new_file)

load = tk.Button(frame1, text='Load config')
load.grid(row = 17, column = 10, sticky='nesw',padx = 10,
                     pady = 0, columnspan=1)
load.bind('<Button-1>', browseFiles)


run = tk.Button(frame1, text='Run Brutus', bg = '#690202', fg = 'white', font='sans 10 bold', command=call_brutus)
run.grid(row = 18, column = 9, sticky='nesw',padx = 10,
                     pady = 4, columnspan=2)
recover('last_used_config.txt')

status0 = tk.Label(frame1, text='Status:')
status0.grid(row=18, column=3, columnspan=1,sticky='ne')

status = tk.Entry(frame1, width=30)
status.grid(row=18, column=4, columnspan=3,sticky='ne',padx=20)

root.mainloop()
