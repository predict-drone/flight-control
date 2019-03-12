# PREDICT (Time-predictable Flight Control)

PREDICT, Time-predictable Flight Control, is a fly-by-wire flight control system based on T-CREST architecture.
The system is connected to X-Plane 11, where you can perform experiemnts on differnet aircrafts with the flight control system. 
The present controller is tuned for Cirrus Vision SF50 aircraft in X-Plane 11. 

One can also implement the flight control on a differnt embedded hardware (like arduino).

The project is open-sourced and all the source codes can be download from this git (except the flight simulator)

![alt text](https://github.com/AAUAutomationAndControl/PREDICT/blob/master/Screenat18_54%20(2).png)

Requirements - 

1. X-Plane 11: Flight-Simulator (https://www.x-plane.com/)
2. X-interface: Interfacing tool (XPC.exe, from the repo) 
3. DE2-115 dev-board with T-CREST (https://github.com/t-crest/patmos) OR other embedded platform.

The interfacing b/w the embedded flight computer and the simulator via x-interface is explained in the following block diagram-

![alt text](https://github.com/AAUAutomationAndControl/PREDICT/blob/master/git_figs.png)



## Setting up X-Interface 

X-Interface is a tool written in C# to interface between X-Plane and embedded flight control computer i.e. T-CREST in this case. However, this can be used to communicate with any embedded hardware (like Arduino) with serial interface. 

```
X-Plane <-------- UDP ---------> X-Interface  <-------- Serial --------->  Embedded Flight Computer (T-CREST)  

```
The default configuration assumes that X-Plane server is at localhost (IP 127.0.0.1) download port: 49001 upload port: 49000. 

![alt text](https://github.com/shibarchi/TiFCon/blob/master/X_plane_settings.PNG)

The default serial port is COM5 (115200 bits/sec, 8-bit word, 1-stop bit, no parity) 
The same can be configured as shown below. Device Manager -> Ports -> USB Serial Port(COM X) -> (right click) -> Properties -> Port settings -> Advanced...

![alt text](https://github.com/shibarchi/TiFCon/blob/master/COM_port.PNG)


## Understanding X-Plane data

X-plane data output can be turned on by opening settings window -> Data Output -> selecting data from the list (Show in Cockpit/ Network via UDP).

For e.g. we have selected the following -

```
Index    Data
  3     Speeds
 17     Pitch...
 26   Throttle(a)
```

The same can be seen as follows in the cockpit.

![alt text](https://github.com/shibarchi/TiFCon/blob/master/cockpit_data.PNG)

The output data format is very interesting and needs understanging of use the data.
We received a string as displayed below - 

![alt text](https://github.com/shibarchi/TiFCon/blob/master/dataread.PNG)

The first 5 elements are ASCII codes (68,65,84,65,42) that represents 'DATA*'.
The next digit is the intex number of the data followed by 3 zeros. i.e. 3,0,0,0 for Speeds.

Next 28 elements represents readings, each feature/ reading is represented by 4 bytes. 
for e.g. the next four elements (106,34,84,67) represents Indicated air speed in knots.

The next 24 elements (,89,12,84,67,67,149,107,67,131,120,106,67,0,192,121,196,181,30,116,67,85,141,135,67,85,141,135,67,) belongs to Vtrue,Vtrue, blank,Vind,Vtru,Vtrue as shown in the cockpit data out. 

The next element is 17 i.e. the index for Pitch,roll, heading followed by 3 zeros. (60,118,187,64,28,193,81,61,1,248,165,66,184,152,135,66,0,192,121,196,0,192,121,196,0,192,121,196,0,192,121,196,) represents pitch, roll, true heading and magnetic heading followed by 4 blank spaces (0,192,121,196 is one blank space). For runnin the sim, we need data from index - 3, 17, 20, 26.

***BitConverter.ToSingle*** (or Int/Int16/Double) can be used in C# to get numeric data from the 4 bytes; 
i.e. **BitConverter.ToSingle(xdata, 45)** will give pitch value if data array is stored as **xdata**


## Bare minimum Control 

To implement a basic autonomous flight control system, we need to control - elevators, ailerons, rudder and throttle. 
In X-Plane control commands for control surfaces are between +1.0 to -1.0 and 0 to 1.0 for throttle. 

## Interfacing
T-CREST can be configured on Altera DE2-115 FPGA devlopment board. Details about configuring T-CREST can be found at https://github.com/t-crest/patmos.

Connect the USB/TTL serial cable (COM5) to the following pins on GPIO:
```
               | * * |  GND
               | * * |
               | * * |
               | * * |
               | * * |  
 RX2-- 39/AH23 | * * |  40/AG26 -- TX2
               +-----+
                GPIO
```
compile and download the flight controller (fcs.c) by typing - 
(make sure the fcs.c file is present in /home/patmos/t-crest/patmos/c/fcs.c)
```
make comp APP=fcs download
```
## Starting the system

1. Set-up T-crest and connect with the host system running the simulator. 
2. Run X-plane, select Cirrus Vision SF50 aircraft and start the simulation.
3. Run X-interface from -> xdirectory-> XPC.exe, and press any key. 
4. Compile and download controller in T-CREST. 

## Understanding the controller

The controller gets sensor data from the simulator via the serial port -

```
290     unsigned char serial =  read(uart2_ptr);
```
The simulator sends a 4 byte header followed by 13 features in 52 byte i.e. 4 bytes per feature. All the parameters are 
in 32-bit floats. The bytes are concatinated to get the parameter value. While concatination, consider the endian of the operating system you are using.  

On successful reception, the state parameters are estimated with a linear quadratic estimator -  
```
333     Estimator( acc_x,  acc_y,  acc_z,  P, Q);
```
### The Estimator 
The Kalman filter maintains the estimates of the state
```
x_hat(k|k)- estimate of x(k) given measurements z(k), z(k-1), ...
x_hat(k+1|k)- estimate of x(k+1) given measurements z(k), z(k-1), ...
```
and the error covariance matrix of the state estimate

```
P(k|k)- covariance of x(k) given z(k), z(k-1), ...
P(k+1|k)- covariance of x(k+1) given z(k), z(k-1), ...
```
State Estimation - 
when, 

x_hat(k|k), u(k) and P(k|k) and measurement z(k+1),
```
State prediction: x_hat(k+1|k) = A(k) x_hat(k|k) + B(k) u(k)
Measurement prediction: z_hat(k+1|1) = C(k) x_hat(k+1|k)

State estimation update: x_hat(k+1|k+1) = x_hat(k+1|k) + W(k+1){z(k+1) - z_hat(k+1|k)}
```
where, W(k+1) is Kalman gain. Read more at https://en.wikipedia.org/wiki/Kalman_filter

### The Controller
The controller is implemented as cascaded PID controllers;

Low-level longitudina and lateral controllers are closed loop controllers with reference signals from High-level controller and feeback signals from the state estimator. 
Low-level controller directly controls the control surfaces to attain reference attitude as commanded by high-lelvel controller. 
High-level controllers are generally open-loop or closed-loop controller that sends command to simulate pilot-inputs or implement predefined flight-modes. 

The structure of the controller can be better understood from the block-diagram below. 

![alt text](https://github.com/AAUAutomationAndControl/PREDICT/blob/master/git_figs_2.png)

 

Now you are ready to take-off!!!
![alt text](https://github.com/AAUAutomationAndControl/PREDICT/blob/master/Sim%20(2).png)
