import math
import os
import matplotlib.patches as patches
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from tkinter import simpledialog
import tkinter as tk
from tkinter import *
from tkinter import messagebox
from functools import partial
import pymongo
import json
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk


class TrajectorySimulator:
    def __init__(self):
        self.__arcs = [
            "Cplongitude",
            "Cplatitude",
            "Start Latitude",
            "End Latitude",
            "Start Longitude",
            "End Longitude",
            "Centre Radius",
            "Clockwise",
            "AntiClockwise",
        ]
        self.__latitude_values = []
        self.__longitude_values = []
        self.__altitude = []
        self.__checkvars = []
        self.__objects = []
        self.__arc_values = []

        self.__data = []
        self.__file = open("Data_Base.json", "w")

        self.root = Tk()
        self.root.geometry("650x850")

        Label(self.root, text="Trajectory Simulator", font="ar 20 bold").grid(
            row=0, column=1
        )

        self.no_of_areas_label = Label(self.root, text="Number of Areas")
        self.no_of_areas_label.grid(row=1, column=0)

        self.no_of_areas_value = StringVar()
        self.no_of_areas_entry = Entry(
            self.root, textvariable=self.no_of_areas_value)
        self.no_of_areas_entry.grid(row=1, column=1)
        Button(self.root, text="Exit", command=self.__self_destroy).grid(
            row=4, column=1
        )

        self.ok_button = Button(self.root, text="OK", command=self.create_area)
        self.ok_button.grid(row=1, column=2)
        self.client = pymongo.MongoClient("mongodb://localhost:27017/")
        self.db = self.client["trajectory_db"]
        self.store_json_button = Button(self.root, text="STORE JSON")
        self.store_json_button.grid(row=2, column=1)
        self.root.mainloop()

    def __self_destroy(self):
        data = {"areas": self.__data}
        json.dump(data, self.__file)
        self.__file.close()
        self.root.destroy()

    def create_area(self):
        num_areas = int(self.no_of_areas_value.get())

        for area in range(num_areas):
            area_window = Toplevel(self.root)
            area_window.geometry("650x850")
            area_window.title(f"Area {area+1}")

            Label(area_window, text="Trajectory Simulator", font="ar 20 bold").grid(
                row=0, column=1
            )

            no_of_coordinates_label = Label(
                area_window, text="Number of Coordinates")
            no_of_coordinates_label.grid(row=1, column=0)

            # Text area :
            no_of_coordinates_value = StringVar()
            no_of_coordinates_entry = Entry(
                area_window, textvariable=no_of_coordinates_value
            )
            no_of_coordinates_entry.grid(row=1, column=1)
            no_of_coordinates_entry.bind(
                "<Return>",
                partial(
                    self.create_coordinates, area, area_window, no_of_coordinates_value
                ),
            )

            # Ok button:
            ok_button = Button(
                area_window,
                text="OK",
                command=partial(
                    self.create_coordinates, area, area_window, no_of_coordinates_value
                ),
            )
            ok_button.grid(row=1, column=2)

            # Submit button:
            submit_button = Button(
                area_window,
                text="Submit",
                command=partial(self.destroy_window, area, area_window),
            )
            submit_button.grid(row=2, column=0, columnspan=2)

    def create_coordinates(self, area, area_window, no_of_cordinates):
        # print(area,' ',area_window,' ',no_of_cordinates.get())
        num_coordinates = int(no_of_cordinates.get())

        Label(area_window, text="Trajectory Simulator", font="ar 20 bold").grid(
            row=0, column=1
        )
        Label(area_window, text="Altitude Range From").grid(row=4, column=0)

        Label(area_window, text="Altitude Range From").grid(row=4, column=0)
        altitudefrom_value = StringVar()
        Entry(area_window, textvariable=altitudefrom_value).grid(row=4, column=1)

        Label(area_window, text="Altitude Range To").grid(row=5, column=0)
        altitudeto_value = StringVar()
        Entry(area_window, textvariable=altitudeto_value).grid(row=5, column=1)

        """ Assigning Variables for future use """
        self.__altitude.append((altitudefrom_value, altitudeto_value))

        latitude_values = []
        longitude_values = []
        checkvars = []
        arc_parameter = []
        for i in range(num_coordinates):
            pos = i * 2 + 7

            # coordinate label
            coordinate = Label(area_window, text=f"Coordinate {i + 1}")
            coordinate.grid(row=pos, column=0)

            # creating checkbox
            checkvar = IntVar()
            checkbtn = Checkbutton(
                area_window,
                text="Arc",
                variable=checkvar,
                command=partial(self.update_parameters,
                                checkvar, i + 1, arc_parameter),
            )
            checkbtn.grid(row=pos, column=1)
            checkvars.append(checkvar)

            # latitude and longitude
            latitude_label = Label(area_window, text="Latitude")
            latitude_label.grid(row=pos + 1, column=0, padx=10)
            latitude_value = StringVar()
            latitude_entry = Entry(area_window, textvariable=latitude_value)
            latitude_entry.grid(row=pos + 1, column=1)

            latitude_values.append(latitude_value)

            longitude_label = Label(area_window, text="Longitude")
            longitude_label.grid(row=pos + 1, column=2, padx=10)
            longitude_value = StringVar()
            longitude_entry = Entry(area_window, textvariable=longitude_value)
            longitude_entry.grid(row=pos + 1, column=3)

            longitude_values.append(longitude_value)

        self.__latitude_values.insert(area, latitude_values)
        self.__longitude_values.insert(area, longitude_values)
        self.__checkvars.insert(area, checkvars)
        self.__arc_values.insert(area, arc_parameter)

    def __crt_file(self, area):
        data = {}
        data["altitude_from"] = self.__altitude[area][0].get()
        data["altitude_to"] = self.__altitude[area][1].get()

        n = len(self.__latitude_values[area])
        data["num_coords"] = n

        k = 0
        cord = []
        for i in range(n):
            dic = {}
            if self.__checkvars[area][i].get() == 1:
                dic["arc"] = "1"
                arc = self.__arc_values[area][k]
                k += 1
                for j in range(9):
                    dic[self.__arcs[j]] = arc[j].get()
            else:
                dic["arc"] = "0"
                dic["latitude"] = self.__latitude_values[area][i].get()
                dic["longitude"] = self.__longitude_values[area][i].get()
            cord.append(dic)
            data["coordinates"] = cord
        print(data)
        self.__data.append(data)
        print(self.__data)

    def destroy_window(self, area, area_window):
        # print("Destroy area : ",area+1)
        # print(f"object{area+1}\n")
        # print(f"altitude from : {self.__altitude[area][0].get()} to : {self.__altitude[area][1].get()} ")
        # for i,each in enumerate(self.__latitude_values[area]):
        #    print(f"latitude{i+1} : {each.get()}")
        # for i,each in enumerate(self.__longitude_values[area]):
        #     print(f"longitude{i+1} : {each.get()}")
        # for j,each in enumerate(self.__arc_values[area]):
        #     print(f"cordinate{j+1} :")
        #     for i,l in enumerate(each):
        #         print(f"element{i+1} : {l.get()}")
        #    print()
        self.__crt_file(area)
        area_window.destroy()

    def update_parameters(self, checkvar, pos, arc_parameter):
        if checkvar.get() == 1:
            # if the arc checkbox is clicked
            self.open_arc_parameters(pos, pos - 1, arc_parameter)
        else:
            # if the checkbox is unchecked
            pass

    def open_arc_parameters(self, coordinate_num, i, arc_parameter):
        arc_window = Toplevel(self.root)
        arc_window.title(f"Arc Coordinate {coordinate_num}")
        arc_window.geometry("300x200")
        parameters = []
        # for all parameters taking input of arc
        Label(arc_window, text="Cplongitude").grid(row=0, column=0)
        cplongitude_value = StringVar()
        Entry(arc_window, textvariable=cplongitude_value).grid(row=0, column=1)
        parameters.append(cplongitude_value)

        Label(arc_window, text="Cplatitude").grid(row=1, column=0)
        cplatitude_value = StringVar()
        Entry(arc_window, textvariable=cplatitude_value).grid(row=1, column=1)
        parameters.append(cplatitude_value)

        Label(arc_window, text="Start Latitude").grid(row=2, column=0)
        startlatitude_value = StringVar()
        Entry(arc_window, textvariable=startlatitude_value).grid(row=2, column=1)
        parameters.append(startlatitude_value)

        Label(arc_window, text="End Latitude").grid(row=3, column=0)
        endlatitude_value = StringVar()
        Entry(arc_window, textvariable=endlatitude_value).grid(row=3, column=1)
        parameters.append(endlatitude_value)

        Label(arc_window, text="Start Longitude").grid(row=4, column=0)
        startlongitude_value = StringVar()
        Entry(arc_window, textvariable=startlongitude_value).grid(row=4, column=1)
        parameters.append(startlongitude_value)

        Label(arc_window, text="End Longitude").grid(row=5, column=0)
        endlongtiude_value = StringVar()
        Entry(arc_window, textvariable=endlongtiude_value).grid(row=5, column=1)
        parameters.append(endlongtiude_value)

        Label(arc_window, text="Centre Radius").grid(row=6, column=0)
        center_rad_value = StringVar()
        Entry(arc_window, textvariable=center_rad_value).grid(row=6, column=1)
        parameters.append(center_rad_value)

        # change here
        clockwise_checkbox_var = BooleanVar(
            value='True')  # Default is clockwise (True)
        anticlockwise_checkbox_var = BooleanVar(
            value='False')  # Default is anticlockwise (False)

        def update_rotation_direction():
            if anticlockwise_checkbox_var.get():
                clockwise_checkbox_var.set('False')
            else:
                clockwise_checkbox_var.set('True')
        Label(arc_window, text="Rotation Direction").grid(row=7, column=0)
        anticlockwise_checkbox = Checkbutton(
            arc_window, text="Anti-clockwise", variable=anticlockwise_checkbox_var, command=update_rotation_direction)
        anticlockwise_checkbox.grid(row=7, column=1)
        parameters.append(clockwise_checkbox_var)
        parameters.append(anticlockwise_checkbox_var)

        submit_button = Button(arc_window, text="Submit",
                               command=arc_window.destroy)
        submit_button.grid(row=8, column=0, columnspan=2)

        arc_parameter.insert(i, parameters)


ob = TrajectorySimulator()

# Create a new figure and set the X and Y labels
#_, ax = plt.subplots()
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Longitude')


# Define the filename for the JSON data
# FILE_NAME = "Data_Base.json"


def do_plotting(FILE_NAME):
    # Check if the file exists at the specified path
    if not os.path.exists(FILE_NAME):
        print(f"File not found at path: {FILE_NAME}")
        return


# Define the filename for the JSON data
FILE_NAME = "Data_Base.json"


def dms_to_decimal(degree_string):
    degree_string = degree_string.strip()
    degrees = float(degree_string[0:degree_string.find("\u00ba")])
    minutes = float(degree_string[degree_string.find(
        "\u00ba") + 1:degree_string.find("'")])
    seconds = float(degree_string[degree_string.find(
        "'") + 1:degree_string.find("\"")])
    decimal_degrees = degrees + (minutes / 60) + (seconds / 3600)
    
    #-----ADDED PORTION---------------------------------------
    direction = degree_string[len(degree_string)-1]
    if direction == 'S' or direction == 'W':
       decimal_degrees = decimal_degrees*-1
    #---------------------------------------------------------
    
    return decimal_degrees


def find_lat_degree(lat):
    return int(lat[0:lat.find("\u00ba")])


def find_long_degree(long):
    return int(long[0:long.find("\u00ba")])


def draw_arc(start_x, start_y, end_x, end_y, center_x, center_y, radius, clock, anti, ax):
    r = float(radius)
    radius_in_degrees = r / 60  # 1 minute = 1 nautical mile = 1/60 degree

    start_angle = math.degrees(math.atan2(
        start_y - center_y, start_x - center_x))
    end_angle = math.degrees(math.atan2(end_y - center_y, end_x - center_x))

    #print("start_angle:", start_angle, "end_angle:", end_angle)
    #print("start_x:", start_x, "start_y:", start_y, "end_x:", end_x, "end_y:", end_y, "center_x:", center_x, "center_y:", center_y)

    #calculated_radius_end = math.sqrt((end_x - center_x)**2 + (end_y - center_y)**2)
    #calculated_radius_start = math.sqrt((start_x - center_x)**2 + (start_y - center_y)**2)

    width = 2 * radius_in_degrees
    height = 2 * radius_in_degrees

    ax.plot([center_x], [center_y], 'ro')

    #print("calculated_radius:", calculated_radius_end, "cal start:", calculated_radius_start, "radius:", radius_in_degrees)

    if not(end_x - start_x > 0 and end_y - start_y > 0):
     #   print("Here")
        t = start_angle
        start_angle = end_angle
        end_angle = t

    arc = patches.Arc((center_x, center_y), width, height, angle=0.0,
                      theta1=start_angle, theta2=end_angle, edgecolor='blue', linewidth=2)
    ax.add_patch(arc)


def draw_arc_2(start_x, start_y, end_x, end_y, ax):
    center_x = (start_x + end_x) / 2.0
    center_y = (start_y + end_y) / 2.0

    # Compute the radius (half of the distance between the endpoints)
    radius = math.sqrt((end_x - start_x)**2 + (end_y - start_y)**2) / 2.0

    # Compute the angle between the x-axis and the line joining the start and end points
    angle = math.degrees(math.atan2(end_y - start_y, end_x - start_x))

    # Draw the semicircular arc from start to end
    arc = patches.Arc((center_x, center_y), 2*radius, 2*radius,
                      angle=angle, theta1=0, theta2=180, edgecolor='black', linewidth=2)
    ax.add_patch(arc)


# List to store the clicked coordinates
#clicked_points = []
#lines_data = []
def on_click(event, ax, clicked_points, lines_data, coordinates_text):
    # Get the coordinates of the click event
    x, y = event.xdata, event.ydata

    # Add the clicked point to the list
    clicked_points.append((x, y))

    # Plot a red point at the clicked location
    ax.plot(x, y, 'ro')

    # Automatically draw lines for the clicked points
    if len(clicked_points) >= 2:
        x_vals, y_vals = zip(*clicked_points)

        # Ask the user to choose between semicircular arc and straight line
        if len(clicked_points) >= 2:
            line_type = simpledialog.askstring(
                "Line Type", "Enter 'arc' for semicircular arc or 'line' for straight line:")
            if line_type == 'arc':
                # Draw a semicircular arc using the two points
                draw_arc_2(x_vals[-2], y_vals[-2], x_vals[-1], y_vals[-1], ax)
                lines_data.append({
                    'type': 'arc',
                    'start': {'latitude': x_vals[-2], 'longitude': y_vals[-2]},
                    'end': {'latitude': x_vals[-1], 'longitude': y_vals[-1]}
                })
            elif line_type == 'line':
                # Draw a straight line connecting the two points
                ax.plot(x_vals[-2:], y_vals[-2:], 'k-')
                lines_data.append({
                    'type': 'line',
                    'start': {'latitude': x_vals[-2], 'longitude': y_vals[-2]},
                    'end': {'latitude': x_vals[-1], 'longitude': y_vals[-1]}
                })

    # Refresh the plot to update the view
    coordinates_text.delete('1.0', tk.END)
    coordinates_text.insert(tk.END, f'({x:.2f}, {y:.2f})')

    plt.draw()


def store_data(clicked_points, lines_data):
    output_format = convert_to_output_format(clicked_points, lines_data)
    # Save the output_format list to a JSON file
    with open('output_format.json', 'w') as f:
        json.dump(output_format, f)
    print("Data stored successfully.")

#store_button = tk.Button(text="Store Data", command=store_data)
# store_button.pack()

# Connect the on_click function to the event handler
#plt.gcf().canvas.mpl_connect('button_press_event', on_click)

north_x=90
north_y=-135
def calculate_angle_between_line_and_semicircle(Cx, Cy, start_x, start_y, x2, y2):
    # Calculate the slope of the line
    slope_line_north = calculate_slope(start_x, start_y, north_x, north_y)
    slope_line=calculate_slope(start_x, start_y, Cx, Cy)
    # Calculate the slope of the tangent line at the point of intersection
    slope_tangent = -1 / slope_line

    # Calculate the angle in radians between the line and the tangent line
    # angle_radians = math.atan(abs(slope_tangent))
    angle_radians = math.atan(abs((slope_tangent - slope_line_north) / (1 + slope_line_north * slope_tangent)))
    # Calculate the angle in degrees
    angle_degrees = math.degrees(angle_radians)

    return angle_degrees
def calculate_slope(x1, y1, x2, y2):
    return (y2 - y1) / (x2 - x1)

def calculate_angle_between_lines(x1, y1, x2, y2, x3, y3):
    # Calculate slopes of the two lines
    m1 = calculate_slope(x1, y1, x2, y2)
    m2 = calculate_slope(x2, y2, x3, y3)

    # Calculate the angle between the two lines
    angle_radians = math.atan(abs((m2 - m1) / (1 + m1 * m2)))
    angle_degrees = math.degrees(angle_radians)

    return angle_degrees

def convert_to_output_format(clicked_points, lines_data):
    output_format = []
    num_points = len(clicked_points)
    output_format.append({
        "SegmentType": "point",
        "StartPoint": {"latitude": clicked_points[0], "longitude": clicked_points[0]},
        "EndPoint": {"latitude": clicked_points[0], "longitude": clicked_points[0]},
        "CourseChange": 0.0,
        "EndAlt": 0.0,
        "StartCourse":0.0,
        "EndCourse": 0.0,
        "Fpa": 0.0,
        "Length": 0.0,
        "TurnCenter": {"latitude": 0.0, "longitude": 0.0},
        "TurnDir": "left To right",
        "TurnRadius": 0.0,
        "windDir_deg": 0,
        "windStr_kts": 0,
        "predictedTAS_kts": 0,
        "ETA": 0,
        "gndSpd_kts": 0.0,
        "predAvail": True,
        "RNP": 0,
        "PredFlightPhase": "UNKNOWN"
    })
    for i in range(len(lines_data)):
        line = lines_data[i]
        if line['type'] == 'line':
            if (i*2 + 1) < num_points:
                start_x, start_y = clicked_points[i*2]
                end_x, end_y = clicked_points[i*2 + 1]
                length = calculate_distance(start_x, start_y, end_x, end_y)
                output_format.append({
                    "SegmentType": "Line",
                    "StartPoint": {"latitude": start_x, "longitude": start_y},
                    "EndPoint": {"latitude": end_x, "longitude": end_y},
                    "CourseChange": 0.0,
                    "EndAlt": 0.0,
                    "StartCourse":calculate_angle_between_lines(north_x,
north_y, start_x, start_y,end_x, end_y),
                    "Fpa": 0.0,
                    "Length": length,
                    "TurnCenter": {"latitude": 0.0, "longitude": 0.0},
                    "TurnDir": "left To right",
                    "TurnRadius": 0.0,
                    "windDir_deg": 0,
                    "windStr_kts": 0,
                    "predictedTAS_kts": 0,
                    "ETA": 0,
                    "gndSpd_kts": 0.0,
                    "predAvail": True,
                    "RNP": 0,
                    "PredFlightPhase": "UNKNOWN"
                })
        elif line['type'] == 'arc':
            if (i*2 + 1) < num_points:
                start_x, start_y = clicked_points[i*2]
                end_x, end_y = clicked_points[i*2 + 1]
                center_x = (start_x + end_x) / 2.0
                center_y = (start_y + end_y) / 2.0
                radius = calculate_distance(
                    center_x, center_y, start_x, start_y)
                output_format.append({
                    "SegmentType": "Arc",
                    "StartPoint": {"latitude": start_x, "longitude": start_y},
                    "EndPoint": {"latitude": end_x, "longitude": end_y},
                    "CourseChange": 0.0,
                    "EndAlt": 0.0,
                    "startCourse": calculate_angle_between_line_and_semicircle(center_x, center_y, start_x, start_y, north_x, north_y),
                    "Fpa": 0.0,
                    "Length": 0.0,
                    "TurnCenter": {"latitude": center_x, "longitude": center_y},
                    "TurnDir": "clockwise",
                    "TurnRadius": radius,
                    "windDir_deg": 0,
                    "windStr_kts": 0,
                    "predictedTAS_kts": 0,
                    "ETA": 0,
                    "gndSpd_kts": 0.0,
                    "predAvail": True,
                    "RNP": 0,
                    "PredFlightPhase": "UNKNOWN"
                })
        output_format.append({
            "SegmentType": "point",
            "StartPoint": {"latitude": clicked_points[-1], "longitude": clicked_points[-1]},
            "EndPoint": {"latitude": clicked_points[-1], "longitude": clicked_points[-1]},
            "CourseChange": 0.0,
            "EndAlt": 0.0,
            "StartCourse":0.0,
            "EndCourse": 0.0,
            "Fpa": 0.0,
            "Length": 0.0,
            "TurnCenter": {"latitude": 0.0, "longitude": 0.0},
            "TurnDir": "left To right",
            "TurnRadius": 0.0,
            "windDir_deg": 0,
            "windStr_kts": 0,
            "predictedTAS_kts": 0,
            "ETA": 0,
            "gndSpd_kts": 0.0,
            "predAvail": True,
            "RNP": 0,
            "PredFlightPhase": "UNKNOWN"
        })
    return output_format


def calculate_distance(start_x, start_y, end_x, end_y):
    return math.sqrt((end_x - start_x)**2 + (end_y - start_y)**2)


def on_zoom(event, ax, canvas):
    ax.set_xlim(event.xdata - 0.1, event.xdata + 0.1)
    ax.set_ylim(event.ydata - 0.1, event.ydata + 0.1)
    canvas.draw()


def zoom_in(ax, canvas):
    current_xlim = ax.get_xlim()
    current_ylim = ax.get_ylim()
    new_xlim = (current_xlim[1] - current_xlim[0]) / 4
    new_ylim = (current_ylim[1] - current_ylim[0]) / 4
    ax.set_xlim(current_xlim[0] + new_xlim, current_xlim[1] - new_xlim)
    ax.set_ylim(current_ylim[0] + new_ylim, current_ylim[1] - new_ylim)
    canvas.draw()


def zoom_out(ax, canvas):
    current_xlim = ax.get_xlim()
    current_ylim = ax.get_ylim()
    new_xlim = (current_xlim[1] - current_xlim[0])
    new_ylim = (current_ylim[1] - current_ylim[0])
    ax.set_xlim(current_xlim[0] - new_xlim, current_xlim[1] + new_xlim)
    ax.set_ylim(current_ylim[0] - new_ylim, current_ylim[1] + new_ylim)
    canvas.draw()

# Your original functions


def draw_area(area, ax):

    n = area['num_coords']
    isPoly = 1

    # Checking whether the area is a straight polygon or not

    for c in area['coordinates']:
        if int(c['arc']) == 1:
            isPoly = 0
            break

    coords_x = []
    coords_y = []
    print("isPoly:", isPoly)

    if isPoly == 1:

        # Drawing straight polygon
        for c in area['coordinates']:
            lat = dms_to_decimal(c['latitude'])
            long = dms_to_decimal(c['longitude'])

            coords_x.append(lat)
            coords_y.append(long)

        coords_x.append(coords_x[0])
        coords_y.append(coords_y[0])

        ax.plot(coords_x, coords_y)
        return

    else:

        # Now it can be either a circle or a section of a circle
        # If it is a circle, then the first and last coordinates will be the same
        if n == 1 and area['coordinates'][0]['arc'] == '1':
            circ = area['coordinates'][0]

            if(circ['Start Latitude'] == circ['End Latitude'] and circ['Start Longitude'] == circ['End Longitude']):
                # It is a circle

                r = float(circ['Centre Radius'])

                x1 = float(dms_to_decimal(circ['Start Latitude']))
                y1 = float(dms_to_decimal(circ['Start Longitude']))

                xc = float(dms_to_decimal(circ['Cplatitude']))
                yc = float(dms_to_decimal(circ['Cplongitude']))

                print("xc:", xc, "yc:", yc, "r:", r)
                circle = plt.Circle((xc, yc), r, fill=False)
                ax.add_patch(circle)

            # else:

                # It is a section of a circle
                # First, draw the straight lines
                #x1 = float(dms_to_decimal(circ['Start Latitude']))
                #y1 = float(dms_to_decimal(circ['Start Longitude']))

                #xc = float(dms_to_decimal(circ['Cplatitude']))
                #yc = float(dms_to_decimal(circ['Cplongitude']))

                #x2 = float(dms_to_decimal(circ['End Latitude']))
                #y2 = float(dms_to_decimal(circ['End Longitude']))

                #ax.plot([x1, x2], [y1, y2])
                #draw_arc(x1, y1, x2, y2, xc, yc, circ['Centre Radius'], circ['Clockwise'], circ['AntiClockwise'], ax)

        for i in range(n):

            c = area['coordinates'][i]
            next_c = area['coordinates'][(i+1) % n]

            # Case 1 : Both normal points
            if(int(c['arc']) == 0) and (int(next_c['arc']) == 0):

                lat = dms_to_decimal(c['latitude'])
                long = dms_to_decimal(c['longitude'])
                next_lat = dms_to_decimal(next_c['latitude'])
                next_long = dms_to_decimal(next_c['longitude'])

                print(lat, long, next_lat, next_long)

                ax.plot([lat, next_lat], [long, next_long])

            # Case 2 : First point is normal, second point is arc
            if(int(c['arc']) == 0) and (int(next_c['arc']) == 1):

                lat = dms_to_decimal(c['latitude'])
                long = dms_to_decimal(c['longitude'])
                next_slat = dms_to_decimal(next_c['Start Latitude'])
                next_slong = dms_to_decimal(next_c['Start Longitude'])

                print(lat, long, next_slat, next_slong)

                ax.plot([lat, next_slat], [long, next_slong])

            # Case 3 : First point is arc, second point is normal
            if(int(c['arc']) == 1) and (int(next_c['arc']) == 0):
                elat = dms_to_decimal(c['End Latitude'])
                elong = dms_to_decimal(c['End Longitude'])
                next_lat = dms_to_decimal(next_c['latitude'])
                next_long = dms_to_decimal(next_c['longitude'])

                print(elat, elong, next_lat, next_long)

                ax.plot([elat, next_lat], [elong, next_long])

            # Case 4 : Both points are arcs
            if(int(c['arc']) == 1) and (int(next_c['arc']) == 1):
                elat = dms_to_decimal(c['End Latitude'])
                elong = dms_to_decimal(c['End Longitude'])
                next_slat = dms_to_decimal(next_c['Start Latitude'])
                next_slong = dms_to_decimal(next_c['Start Longitude'])

                print(elat, elong, next_slat, next_slong)

                plt.plot([elat, next_slat], [elong, next_slong])

            # Finally, draw the arcs
            if int(c['arc']) == 1:
                slat = dms_to_decimal(c['Start Latitude'])
                elat = dms_to_decimal(c['End Latitude'])
                slong = dms_to_decimal(c['Start Longitude'])
                elong = dms_to_decimal(c['End Longitude'])
                cplat = dms_to_decimal(c['Cplatitude'])
                cplong = dms_to_decimal(c['Cplongitude'])

                print(slat, slong, elat, elong, cplat, cplong)

                draw_arc(slat, slong, elat, elong, cplat, cplong,
                         c['Centre Radius'], c['Clockwise'], c['AntiClockwise'], ax)

# Load the JSON data from the file
#f = open('./'+FILE_NAME)
#data = json.load(f)
#area_arr = data['areas']

# Loop through each area in the JSON data and draw it on the plot
# for area in area_arr:
    #draw_area(area, ax)

# Set the aspect ratio to 'equal' and display the grid
# lt.gca().set_aspect('equal')
# plt.grid(True)

# Show the plot
# plt.show()


def do_plotting(FILE_NAME):

    # Create a new figure and set the X and Y labels
    fig, ax = plt.subplots()
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Longitude')

    root = tk.Tk()
    root.title("Trajectory Plot")

    coordinates_text = tk.Text(root, height=1, width=15, wrap=tk.NONE)
    coordinates_text.pack(side=tk.BOTTOM)

    # Load the JSON data from the file
    f = open('./'+FILE_NAME)
    data = json.load(f)
    area_arr = data['areas']

    # Loop through each area in the JSON data and draw it on the plot
    for area in area_arr:
        draw_area(area, ax)

    # List to store the clicked coordinates
    clicked_points = []
    lines_data = []

    # lambda is used to pass arguments to the function

    # Connect the on_click function to the event handler
    plt.gcf().canvas.mpl_connect('button_press_event', lambda event: on_click(
        event, ax, clicked_points, lines_data, coordinates_text))

    # Set the aspect ratio to 'equal' and display the grid
    plt.gca().set_aspect('equal')
    plt.grid(True)

    # Show the plot
    plt.ion()

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    hscrollbar = tk.Scrollbar(root, orient=tk.HORIZONTAL,
                              command=canvas.get_tk_widget().xview)
    hscrollbar.pack(side=tk.BOTTOM, fill=tk.X)
    canvas.get_tk_widget().config(xscrollcommand=hscrollbar.set)

    # Create vertical scrollbar
    vscrollbar = tk.Scrollbar(root, orient=tk.VERTICAL,
                              command=canvas.get_tk_widget().yview)
    vscrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    canvas.get_tk_widget().config(yscrollcommand=vscrollbar.set)

    plot_frame = tk.Frame(canvas.get_tk_widget())
    canvas.get_tk_widget().create_window((0, 0), window=plot_frame, anchor=tk.NW)

    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    canvas.mpl_connect('button_press_event',
                       lambda event: on_zoom(event, ax, canvas))

    store_button = tk.Button(
        text="Store Data", command=lambda: store_data(clicked_points, lines_data))
    store_button.pack(side=tk.TOP)

    # Create "Zoom In" button
    zoom_in_button = tk.Button(
        root, text="Zoom In", command=lambda: zoom_in(ax, canvas))
    zoom_in_button.pack(side=tk.LEFT)

    # Create "Zoom Out" button
    zoom_out_button = tk.Button(
        root, text="Zoom Out", command=lambda: zoom_out(ax, canvas))
    zoom_out_button.pack(side=tk.LEFT)

    toolbar = NavigationToolbar2Tk(canvas, plot_frame)
    toolbar.update()
    toolbar.pack(side=tk.TOP, fill=tk.X)

    root.mainloop()


do_plotting(FILE_NAME)
# # to be updated
# north_x=0
# north_y=0
# def calculate_angle_between_line_and_semicircle(Cx, Cy, start_x, start_y, x2, y2):
#     # Calculate the slope of the line
#     slope_line_north = calculate_slope(start_x, start_y, north_x, north_y)
#     slope_line=calculate_slope(start_x, start_y, Cx, Cy)
#     # Calculate the slope of the tangent line at the point of intersection
#     slope_tangent = -1 / slope_line

#     # Calculate the angle in radians between the line and the tangent line
#     # angle_radians = math.atan(abs(slope_tangent))
#     angle_radians = math.atan(abs((slope_tangent - slope_line_north) / (1 + slope_line_north * slope_tangent)))
#     # Calculate the angle in degrees
#     angle_degrees = math.degrees(angle_radians)

#     return angle_degrees
# def calculate_slope(x1, y1, x2, y2):
#     return (y2 - y1) / (x2 - x1)

# def calculate_angle_between_lines(x1, y1, x2, y2, x3, y3):
#     # Calculate slopes of the two lines
#     m1 = calculate_slope(x1, y1, x2, y2)
#     m2 = calculate_slope(x2, y2, x3, y3)

#     # Calculate the angle between the two lines
#     angle_radians = math.atan(abs((m2 - m1) / (1 + m1 * m2)))
#     angle_degrees = math.degrees(angle_radians)

#     return angle_degrees