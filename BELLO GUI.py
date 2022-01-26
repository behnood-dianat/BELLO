import PySimpleGUI as sg
import Angle_sorter
import BELLO_main
import Radial_Pair_Distribution_Function
import Coordination_Heatmap


def true_call_string(values):
    temp=[]
    for x in values:
        if x=='' or x==' ':
            temp.append(False)
        else:
            for i in str(x):
                if i !=' ' and i !='':
                    temp.append(True)
                else:
                    print('false is here',i)
                    temp.append(False)
    if all(temp) and len(temp) > 1:
        return (True)

def true_call(values):
    temp=[]
    for x in values:
        if x=='' or x==' ':
            temp.append(False)
        else:
            for i in x:
                if i in ('0123456789.') and i !='':
                    temp.append(True)
                else:
                    temp.append(False)
    if all(temp) and len(temp) > 1:
        return(True)

def window_rdf():
    RDF_value= True
    layout_rdf=[
            [sg.Push(),sg.Text('Maximum r'),sg.InputText(default_text='',key='r_max',size=(20,10))],
            [sg.Push(),sg.Text('delta r'),sg.InputText(default_text='',key='delta_r',size=(20,10))],
            [sg.Push(), sg.Text('First element'), sg.InputText(default_text='',key='First_element', size=(20, 10))],
            [sg.Push(), sg.Text('Second element'), sg.InputText(default_text='',key='Second_element', size=(20, 10))],
            [sg.Push(),sg.Button('Submit',key="submit")],
            ]
    window_rdf= sg.Window('RDF input',layout_rdf,keep_on_top=True)
    while True:
        event, values = window_rdf.read()
        if event=='submit' and values!='' and true_call(values['r_max']) and true_call(values['delta_r'])\
                and true_call_string(values['First_element']) and true_call_string(values['Second_element']):
            print('Values submitted')
            window_rdf.close()
            return(values['delta_r'],values['r_max'],values['First_element'],values['Second_element'],RDF_value)

        elif event==sg.WINDOW_CLOSED:
            RDF_value = False
            return([RDF_value])

def window_angle():
    angle_value= True
    def make_layout(number):
        if number==2:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='2')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==3:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='3')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==4:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='4')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Text('Element_4'),sg.Push(),sg.InputText(default_text='',key='Element_four',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==5:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='5')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Text('Element_4'),sg.Push(),sg.InputText(default_text='',key='Element_four',size=(5,10))],
                    [sg.Text('Element_5'),sg.Push(),sg.InputText(default_text='',key='Element_fifth',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==6:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='6')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Text('Element_4'),sg.Push(),sg.InputText(default_text='',key='Element_four',size=(5,10))],
                    [sg.Text('Element_5'),sg.Push(),sg.InputText(default_text='',key='Element_fifth',size=(5,10))],
                    [sg.Text('Element_6'),sg.Push(),sg.InputText(default_text='',key='Element_sixth',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        return(layout)
    layout=[
            [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'),default_value='', size=(5, 1),enable_events=True)],
            [sg.Push(),sg.Button('Submit',key="submit",enable_events=True)],
            ]
    angle_window= sg.Window('Angle sorting input',layout,keep_on_top=True)
    while True:
        event, values = angle_window.read()
        if event=='submit' and true_call_string(list(values.values())[1:len(values)]):
            angle_window.close()
            a=list(values.values())
            print('Values submitted')
            return(angle_value,a[1:len(a)])
        if event==sg.WINDOW_CLOSED:
            angle_value = False
            return(angle_value,list({'x': 'x'}.values()))
        if values[0]=='2':
            angle_window.close()
            angle_window = sg.Window('Angle sorting input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='3':
            angle_window.close()
            angle_window = sg.Window('Angle sorting input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='4':
            angle_window.close()
            angle_window = sg.Window('Angle sorting input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='5':
            angle_window.close()
            angle_window = sg.Window('Angle sorting input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='6':
            angle_window.close()
            angle_window = sg.Window('Angle sorting input',make_layout(int(values[0])),keep_on_top=True)

def window_coord():
    coord_value= True
    def make_layout(number):
        if number==2:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='2')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==3:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='3')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==4:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='4')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Text('Element_4'),sg.Push(),sg.InputText(default_text='',key='Element_four',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==5:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='5')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Text('Element_4'),sg.Push(),sg.InputText(default_text='',key='Element_four',size=(5,10))],
                    [sg.Text('Element_5'),sg.Push(),sg.InputText(default_text='',key='Element_fifth',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        elif number==6:
            layout=[
                    [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'), size=(5, 1),enable_events=True,default_value='6')],
                    [sg.Text('Element_1'),sg.Push(),sg.InputText(key='Element_one',size=(5,10),default_text='Ge')],
                    [sg.Text('Element_2'),sg.Push(),sg.InputText(default_text='',key='Element_two',size=(5,10))],
                    [sg.Text('Element_3'),sg.Push(),sg.InputText(default_text='',key='Element_three',size=(5,10))],
                    [sg.Text('Element_4'),sg.Push(),sg.InputText(default_text='',key='Element_four',size=(5,10))],
                    [sg.Text('Element_5'),sg.Push(),sg.InputText(default_text='',key='Element_fifth',size=(5,10))],
                    [sg.Text('Element_6'),sg.Push(),sg.InputText(default_text='',key='Element_sixth',size=(5,10))],
                    [sg.Push(),sg.Button('Submit',key="submit")],
                    ]
        return(layout)
    layout=[
            [sg.Text('Number of chemical elements'),sg.Push(),sg.Combo(('2', '3','4','5','6'),default_value='', size=(5, 1),enable_events=True)],
            [sg.Push(),sg.Button('Submit',key="submit",enable_events=True)],
            ]
    coord_window= sg.Window('Coordination number chemical specie input',layout,keep_on_top=True)
    while True:
        event, values = coord_window.read()
        if event=='submit' and true_call_string(list(values.values())[1:len(values)]):
            coord_window.close()
            a=list(values.values())
            print('Values submitted')
            return(coord_value,a[1:len(a)])
        if event==sg.WINDOW_CLOSED:
            coord_value = False
            return(coord_value,list({'x': 'x'}.values()))
        if values[0]=='2':
            coord_window.close()
            coord_window = sg.Window('Coordination number input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='3':
            coord_window.close()
            coord_window = sg.Window('Coordination number input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='4':
            coord_window.close()
            coord_window = sg.Window('Coordination number input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='5':
            coord_window.close()
            coord_window = sg.Window('Coordination number input',make_layout(int(values[0])),keep_on_top=True)
        if values[0]=='6':
            coord_window.close()
            coord_window = sg.Window('Coordination number input',make_layout(int(values[0])),keep_on_top=True)

sg.theme("darkamber")
layout_main = [
    [   [sg.Image('cover.png')],
        [sg.Text('Insert xyz file')],
        [sg.Input(key='url',expand_x=True,default_text = ''),
         sg.FileBrowse(target=(2, 0),file_types=(("XYZ Files", "*.xyz"), ("ALL Files", "*.*")))]],
        [sg.Text('Parameters'),sg.HorizontalSeparator()],
        [sg.Text('Inter-atomic distance threshold'),
         sg.Checkbox('Automatic',enable_events=True, default=False, key="automatic"),sg.Push(),
         sg.InputText(key="trh",size=(10,10),default_text='3.0')],
        [sg.Text('Inter-atomic distance tolerance'),sg.Push(),sg.InputText(key="tlr",size=(10,10),default_text='0.5')],
        [sg.Text('Separate Angle distributions (memory efficient)'),sg.Push(),
         sg.Checkbox('',enable_events=True, default=False, key="sep_ang")],
        [sg.Text(' ')],
        [sg.Text('Unitcell dimensions (Ang)'),sg.Push(),sg.Text('X'),sg.Input(key='x',size=(20,10),default_text='20.00000')],
        [sg.Push(),sg.Text('Y'),sg.Input(key='y',size=(20,10),default_text='20.00000')],
        [sg.Push(),sg.Text('Z'),sg.Input(key='z',size=(20,10),default_text='20.00000')],
        [sg.Text('Calculations'),sg.HorizontalSeparator()],
        [sg.Radio('BELLO','radios',               enable_events=True, default=False, key="BELLO"),
         sg.Radio('Coordination Heatmap','radios',enable_events=True, default=False, key="COORD"),
         sg.Radio('Angle sorter','radios',        enable_events=True, default=False, key="ANGLE"),
         sg.Radio('RDF','radios',                 enable_events=True, default=False, key="RDF")],
        [sg.Push(),sg.Button('Calculate',                                            key="calculate")],
]

key_values=['BELLO','COORD','ANGLE','RDF']
window = sg.Window('BELLO Input', layout_main)

while True:
    event, values = window.read()
    if event==sg.WINDOW_CLOSED:
        break
    if values["automatic"]==True and event=='automatic':
        window['trh'].update(disabled=True)
        window['tlr'].update(disabled=True)
    if values["automatic"]==False and event=='automatic':
        window['trh'].update(disabled=False)
        window['tlr'].update(disabled=False)

    if event=='calculate' \
            and values['BELLO']==True \
            and values['url'] != '' \
            and true_call([values['x'],values['y'],values['z']]) \
            and values["automatic"]==True:
        print("Running BELLO (Automatic threshold)")
        BELLO_main.BELLO(values['url'],values['x'],values['y'],values['z'],values['automatic'],2,1,values['sep_ang'])

    if event=='calculate' \
            and values['BELLO']==True \
            and values['url'] != '' \
            and true_call([values['x'],values['y'],values['z']]) \
            and values["automatic"]==False:
        print("Running BELLO")
        BELLO_main.BELLO(values['url'],values['x'],values['y'],values['z'],values['automatic'],values['trh'],values['tlr'],values['sep_ang'])

    if (event=='calculate'
        and values['url']==''
        and (values['BELLO'] == True or values['RDF'] == True)):
        print("Need XYZ file")

    if (event=='calculate'
        and (values['BELLO']==True or values['RDF']==True)
        and true_call([values['x'],values['y'],values['z']])==False):
        print("Need float Cell dimensions")

    if event=='RDF' and values['RDF']==True:
        RDF_inputs = window_rdf()
        window['RDF'].update(RDF_inputs[-1])
    elif event=='RDF' and values['RDF']==False:
        print('RDF: off')

    if event == 'COORD' and values['COORD']==True:
        values['COORD'], coord_elements = window_coord()
        window['COORD'].update(values['COORD'])

    elif event=='COORD' and values['COORD']==False:
        print('Cordination Heatmap: off')

    if event=='ANGLE' and values['ANGLE']==True:
        values['ANGLE'], angle_elements = window_angle()
        window['ANGLE'].update(values['ANGLE'])

    elif event=='ANGLE' and values['ANGLE']==False:
        print('Angle sorter: off')

    if event=='calculate' \
            and values['ANGLE']==True:
        Angle_sorter.sorter(angle_elements)

    if event=='calculate' \
            and values['COORD']==True:
        Coordination_Heatmap.coordination_heatmap(coord_elements)

    if event=='calculate' \
        and values['url']!='' \
        and values['RDF']==True \
        and true_call([values['x'],values['y'],values['z']]):
        if len(RDF_inputs) > 2:
            Radial_Pair_Distribution_Function.RDF(values['url'],values['x'],values['y'],values['z'],
                                                  RDF_inputs[0],RDF_inputs[1],RDF_inputs[2],RDF_inputs[3])