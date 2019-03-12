//
//   Copyright 2018 Aalborg University, Denamrk. 
//   All rights reserved.
//   Author: Shibarchi Majumder (sm@es.aau.dk)
//		     Oktay Baris (okba@dtu.dk)
//   
//   
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are met:
//
//     1. Redistributions of source code must retain the above copyright notice,
//        author, this list of conditions and the following disclaimer.
//
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, author, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ``AS IS'' AND ANY EXPRESS
//   OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
//   NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
//   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
//   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   The views and conclusions contained in the software and documentation are
//   those of the authors and should not be interpreted as representing official
//   policies, either expressed or implied, of the copyright holder.
//
//
const int NOC_MASTER = 0;
#include <machine/patmos.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "libcorethread/corethread.h"

#define AUTO_TAKEOFF  1
#define ALT_HOLD      2
#define AUTO_THROTTLE 3
#define CLIMB         4


typedef unsigned char BYTE;

_iodev_ptr_t uart1_ptr = (_iodev_ptr_t) 0xF0080000;
_iodev_ptr_t uart2_ptr = (_iodev_ptr_t) 0xF00E0000;
_iodev_ptr_t uart3_ptr = (_iodev_ptr_t) 0xF00F0000;

BYTE out_elev_H, out_elev_L, out_ail_H, out_ail_L, out_rud_H, out_rud_L, out_thr_H, out_thr_L;

int mode = AUTO_TAKEOFF;


float pitch = 0;
float pitch2 = 0;
float roll = 0;
float heading = 0;
const float dt = 0.05;
float longitude, latitude;

// Data Fusion parameters 

_UNCACHED volatile float acc_x = 0;
_UNCACHED volatile float acc_y = 0;
_UNCACHED volatile float acc_z = 0;
_UNCACHED volatile float P = 0;
_UNCACHED volatile float Q = 0;

_UNCACHED volatile float ias_sen = 0;
_UNCACHED volatile float p_ias_sen = 0;
_UNCACHED volatile float alt_ft_msl = 0;
_UNCACHED volatile float alt_ft_agl = 0; 

_UNCACHED volatile float theta_cal = 0.0;
_UNCACHED volatile float phi_cal = 0.0; 

float P0_1_1 = 0;
float P0_2_1 = 0;
float P0_1_2 = 0;
float P0_2_2 = 0;
float P1_1_1 = 0;
float P1_2_1 = 0;
float P1_1_2 = 0;
float P1_2_2 = 0;

const float gyro_noise = 0.001;
const float gyro_bias = 0.03;
const float Sensor_accuracy = 0.03;

float Kalman_gain_0_1 = 0;
float Kalman_gain_0_2 = 0;
float Kalman_gain_1_1 = 0;
float Kalman_gain_1_2 = 0;

float bias_1 = 0;
float bias_2 = 0;

float c99max (float a, float lowest){
	if (a< lowest){
		a = lowest;
	}
	return a; 
}



//Thread 1
void Estimator1(void *arg) {


		acc_x = (acc_x - (ias_sen - p_ias_sen)/dt);
		float theta_acc = (float)(atan2(acc_x, sqrt(acc_z*acc_z + acc_y*acc_y))*57.324);

		float temp1 = P - bias_1;
		theta_cal += dt * temp1;

		P0_1_1 += dt * (dt* P0_2_2 - P0_1_2 - P0_2_1 + gyro_noise);
		P0_1_2 -= dt * P0_2_2;
		P0_2_1 -= dt * P0_2_2;
		P0_2_2 += gyro_bias * dt;

		float temp3 = c99max(P0_1_1 + Sensor_accuracy, 0.00001);

		//calculate Kalman Gain
		Kalman_gain_0_1 = c99max(P0_1_1 / temp3, .000001);
		Kalman_gain_0_2 = c99max(P0_1_2 / temp3, .000001);


		// calculate system state
		theta_cal += Kalman_gain_0_1 * (theta_acc - theta_cal);
		bias_1 += Kalman_gain_0_2 * (theta_acc - theta_cal);	

		//Calculate Aposteriori Covarience
		P0_1_1 -= Kalman_gain_0_1 * P0_1_1;
	    P0_1_2 -= Kalman_gain_0_1 * P0_1_2;
	    P0_2_1 -= Kalman_gain_0_2 * P0_2_1;
	    P0_2_2 -= Kalman_gain_0_2 * P0_2_2;


}



// Thread 2
void Estimator2(void *arg) {

		p_ias_sen = ias_sen;
		float phi_acc = (float)(atan2(acc_y,acc_z)*57.324);

		float temp2 = Q- bias_2;
		phi_cal += dt * temp2;

		P1_1_1 += dt * (dt* P1_2_2 - P1_1_2 - P1_2_1 + gyro_noise);
		P1_1_2 -= dt * P1_2_2;
		P1_2_1 -= dt * P1_2_2;
		P1_2_2 += gyro_bias * dt;

		float temp4 = c99max(P1_1_1 + Sensor_accuracy, 0.00001);

		//calculate Kalman Gain
		Kalman_gain_1_1 = c99max(P1_1_1 / temp4, .000001);
		Kalman_gain_1_2 = c99max(P1_1_2 / temp4, .000001);


		// calculate system state
		phi_cal += Kalman_gain_1_1 * (phi_acc - phi_cal);	
		bias_2 += Kalman_gain_1_2 * (phi_acc - phi_cal);

		//Calculate Aposteriori Covarience
		P1_1_1 -= Kalman_gain_1_1 * P1_1_1;
	    P1_1_2 -= Kalman_gain_1_1 * P1_1_2;
	    P1_2_1 -= Kalman_gain_1_2 * P1_2_1;
	    P1_2_2 -= Kalman_gain_1_2 * P1_2_2;

}




void Estimator(void *arg){

    acc_x = (acc_x - (ias_sen - p_ias_sen)/dt);
	p_ias_sen = ias_sen;

	float theta_acc = (float)(atan2(acc_x, sqrt(acc_z*acc_z + acc_y*acc_y))*57.324);
	float phi_acc = (float)(atan2(acc_y,acc_z)*57.324);
	
	float temp1 = P - bias_1;
    theta_cal += dt * temp1;

	float temp2 = Q- bias_2;
	phi_cal += dt * temp2;

	P0_1_1 += dt * (dt* P0_2_2 - P0_1_2 - P0_2_1 + gyro_noise);
	P0_1_2 -= dt * P0_2_2;
	P0_2_1 -= dt * P0_2_2;
	P0_2_2 += gyro_bias * dt;

	P1_1_1 += dt * (dt* P1_2_2 - P1_1_2 - P1_2_1 + gyro_noise);
	P1_1_2 -= dt * P1_2_2;
	P1_2_1 -= dt * P1_2_2;
	P1_2_2 += gyro_bias * dt;

	float temp3 = c99max(P0_1_1 + Sensor_accuracy, 0.00001);
	float temp4 = c99max(P1_1_1 + Sensor_accuracy, 0.00001);

	Kalman_gain_0_1 = c99max(P0_1_1 / temp3, .000001);
	Kalman_gain_0_2 = c99max(P0_1_2 / temp3, .000001);
	Kalman_gain_1_1 = c99max(P1_1_1 / temp4, .000001);
	Kalman_gain_1_2 = c99max(P1_1_2 / temp4, .000001);

	theta_cal += Kalman_gain_0_1 * (theta_acc - theta_cal);
	phi_cal += Kalman_gain_1_1 * (phi_acc - phi_cal);
	bias_1 += Kalman_gain_0_2 * (theta_acc - theta_cal);	
	bias_2 += Kalman_gain_1_2 * (phi_acc - phi_cal);

	P0_1_1 -= Kalman_gain_0_1 * P0_1_1;
    P0_1_2 -= Kalman_gain_0_1 * P0_1_2;
    P0_2_1 -= Kalman_gain_0_2 * P0_2_1;
    P0_2_2 -= Kalman_gain_0_2 * P0_2_2;
	P1_1_1 -= Kalman_gain_1_1 * P1_1_1;
    P1_1_2 -= Kalman_gain_1_1 * P1_1_2;
    P1_2_1 -= Kalman_gain_1_2 * P1_2_1;
    P1_2_2 -= Kalman_gain_1_2 * P1_2_2;

}


int receiving = 0;
int bit_1 = 0;
int bit_2 = 0;
int bit_3 = 0;

int elev, ail, rudder, throttle, n_wheel; 

int all_received = 0;

int total_bytes = 52;

int i = 0;
unsigned char raw_bytes_from_sim[53];//original value was 52

float p_control_accu_err = 0;
float p_control_d_err = 0;
float cruise_control_accu_err = 0;
float cruise_control_d_err = 0; 


float pitch_control_prop_gain = 2;
float pitch_control_diff_gain = 2;
float pitch_control_int_gain = 0.001;

float cruise_control_prop_gain = 0.020;
float cruise_control_diff_gain = 0.02;
float cruise_control_int_gain = 0.0001;



void longitudinal_control(int d_pitch) {
    rudder = 0;
	float err = d_pitch - pitch;
	p_control_accu_err += err;
	float der = err - p_control_d_err;
	p_control_d_err = err;
	elev = (int) (err * pitch_control_prop_gain + (p_control_accu_err *pitch_control_int_gain) + der * pitch_control_diff_gain);
	//printf("err: %f, derr: %f, ierr: %f, elev: %d \n", err, der, p_control_accu_err, elev );

}


float min_max(float max, float min, float val){
	if (val > max){
		val = max;
	}
	else if(val < min){
		val = min;
	}
	
	return val;
}


void cruise_control(int altitude) {
    rudder = 0;
	throttle = 60;
	float err = altitude - alt_ft_msl;
	cruise_control_accu_err += err;
	float der = err - cruise_control_d_err;

	longitudinal_control(min_max(15,-15,(err * cruise_control_prop_gain )));
	
}

void delay(int a) {

	for (int i = 0; i < a * 10000; i++) {
	
	}
}


void inline write(_iodev_ptr_t uart_base_ptr, unsigned char c) {	
  while ((*uart_base_ptr & 0x01) == 0);
  *(uart_base_ptr+1) = c;
}

char inline read(_iodev_ptr_t uart_base_ptr) {		
  while ((*uart_base_ptr & 0x02) == 0);
  return *(uart_base_ptr+1);
}

int safe_byte(float a) {
	int temp = (int) a; 
	temp = temp + 100;

	if (temp > 200) {
		temp = 200;
	}
	else if (temp < 0) {
		temp = 0;
	}

	return temp;
}

float roll_p = 1;
float roll_d = 30;
float roll_i = 0.0001;
float r_control_accu_err = 0.0;
float r_control_d_err = 0.0;



void lateral_control(float d_roll) {

    float err = d_roll - roll;
	r_control_accu_err += err;
	float der = err - r_control_d_err;
	r_control_d_err = err;
	ail = (int) (err * roll_p + (r_control_accu_err *roll_i) + der * roll_d);

}



void HDG_control (float d_hdg){
 float err = d_hdg - heading;
   if(err >= 180){
	   lateral_control(min_max(10,-10,-err*1));
   }
   else{
	   lateral_control(min_max(10,-10,err*1));
   }
}

float conv_to_float(unsigned char byte1, unsigned char byte2,unsigned char byte3, unsigned char byte4 ){

	union u_tag {
		unsigned char b[4];
		float fval;
	} u;

	u.b[0] = byte1;
	u.b[1] = byte2;
	u.b[2] = byte3;
	u.b[3] = byte4;

	return u.fval;

}



int main() {
    printf("Patmos Started!!");



	//int id = get_cpuid();
	//int cnt = get_cpucnt();

	for(;;){

		unsigned char serial =  read(uart2_ptr);

		if (receiving == 1) {
			i += 1;
			
		}
		if (serial == 255 && receiving == 0) {
			bit_1 = 1;
		}
		else if (serial == 254 && receiving == 0 && bit_1 == 1){
			bit_2 = 1;
		}
		else if (serial == 253 && receiving == 0 && bit_2 == 1) {
			bit_3 = 1;
		}
		else if (serial == 252 && receiving == 0 && bit_3 == 1) {
			receiving = 1;
			bit_1 = 0;
			bit_2 = 0;
			bit_3 = 0;
		}
		else {
			bit_1 = 0;
			bit_2 = 0;
			bit_3 = 0;
		}

		raw_bytes_from_sim[i] = serial;

		if (i == total_bytes) {
				receiving = 0;
				i = 0;
				bit_1 = 0;
			    bit_2 = 0;
			    bit_3 = 0;
				all_received = 1;
		}


        if(all_received == 1){




				//unsigned char b2[] = { raw_bytes_from_sim[4], raw_bytes_from_sim[3],raw_bytes_from_sim[2],raw_bytes_from_sim[1] };
				//memcpy(&pitch, &b2, sizeof(pitch));

        		//pitch2 = (float)((raw_bytes_from_sim[4] << 24) | (raw_bytes_from_sim[3] << 16) | (raw_bytes_from_sim[2] <<8) | (raw_bytes_from_sim[1])) ;

				pitch = conv_to_float(raw_bytes_from_sim[4], raw_bytes_from_sim[3],raw_bytes_from_sim[2],raw_bytes_from_sim[1]);

				//printf(" Pitch1 : %b \n", pitch);
				//printf(" Pitch2 : %b \n", pitch2);


				//unsigned char c[] = { raw_bytes_from_sim[8], raw_bytes_from_sim[7],raw_bytes_from_sim[6],raw_bytes_from_sim[5] };
				//memcpy(&roll, &c, sizeof(roll));

				//roll = (raw_bytes_from_sim[5] << 24) | (raw_bytes_from_sim[6] << 16) | (raw_bytes_from_sim[7] <<8) | (raw_bytes_from_sim[8]) ;
				roll = conv_to_float(raw_bytes_from_sim[8], raw_bytes_from_sim[7],raw_bytes_from_sim[6],raw_bytes_from_sim[5]);

				//unsigned char d[] = { raw_bytes_from_sim[12], raw_bytes_from_sim[11],raw_bytes_from_sim[10],raw_bytes_from_sim[9] };
				//memcpy(&heading, &d, sizeof(heading));

				//heading = (raw_bytes_from_sim[9] << 24) | (raw_bytes_from_sim[10] << 16) | (raw_bytes_from_sim[11] <<8) | (raw_bytes_from_sim[12]) ;
				heading = conv_to_float(raw_bytes_from_sim[12], raw_bytes_from_sim[11],raw_bytes_from_sim[10],raw_bytes_from_sim[9]);

				//unsigned char e[] = { raw_bytes_from_sim[16], raw_bytes_from_sim[15],raw_bytes_from_sim[14],raw_bytes_from_sim[13] };
				//memcpy(&ias_sen, &e, sizeof(ias_sen));

				//ias_sen = (raw_bytes_from_sim[13] << 24) | (raw_bytes_from_sim[14] << 16) | (raw_bytes_from_sim[15] <<8) | (raw_bytes_from_sim[16]) ;
				ias_sen = conv_to_float(raw_bytes_from_sim[16], raw_bytes_from_sim[15],raw_bytes_from_sim[14],raw_bytes_from_sim[13]);

				//unsigned char f[] = { raw_bytes_from_sim[20], raw_bytes_from_sim[19],raw_bytes_from_sim[18],raw_bytes_from_sim[17] };
				//memcpy(&alt_ft_msl, &f, sizeof(alt_ft_msl));

				//alt_ft_msl = (raw_bytes_from_sim[17] << 24) | (raw_bytes_from_sim[18] << 16) | (raw_bytes_from_sim[19] <<8) | (raw_bytes_from_sim[20]) ;
				alt_ft_msl = conv_to_float(raw_bytes_from_sim[20], raw_bytes_from_sim[19],raw_bytes_from_sim[18],raw_bytes_from_sim[17]);

				//unsigned char g[] = { raw_bytes_from_sim[24], raw_bytes_from_sim[23],raw_bytes_from_sim[22],raw_bytes_from_sim[21] };
				//memcpy(&alt_ft_agl, &g, sizeof(alt_ft_agl));

				//alt_ft_agl = (raw_bytes_from_sim[21] << 24) | (raw_bytes_from_sim[22] << 16) | (raw_bytes_from_sim[23] <<8) | (raw_bytes_from_sim[24]) ;
				alt_ft_agl = conv_to_float(raw_bytes_from_sim[24], raw_bytes_from_sim[23],raw_bytes_from_sim[22],raw_bytes_from_sim[21]);

				//unsigned char h[] = { raw_bytes_from_sim[28], raw_bytes_from_sim[27],raw_bytes_from_sim[26],raw_bytes_from_sim[25] };
				//memcpy(&acc_x, &h, sizeof(acc_x));

				//acc_x = (raw_bytes_from_sim[25] << 24) | (raw_bytes_from_sim[26] << 16) | (raw_bytes_from_sim[27] <<8) | (raw_bytes_from_sim[28]) ;
				acc_x = conv_to_float(raw_bytes_from_sim[28], raw_bytes_from_sim[27],raw_bytes_from_sim[26],raw_bytes_from_sim[25]);

				
				//unsigned char i[] = { raw_bytes_from_sim[32], raw_bytes_from_sim[31],raw_bytes_from_sim[30],raw_bytes_from_sim[29] };
				//memcpy(&acc_y, &i, sizeof(acc_y));

				//acc_y = (raw_bytes_from_sim[29] << 24) | (raw_bytes_from_sim[30] << 16) | (raw_bytes_from_sim[31] <<8) | (raw_bytes_from_sim[32]) ;
				acc_y = conv_to_float(raw_bytes_from_sim[32], raw_bytes_from_sim[31],raw_bytes_from_sim[30],raw_bytes_from_sim[29]);


				//unsigned char j[] = { raw_bytes_from_sim[36], raw_bytes_from_sim[35],raw_bytes_from_sim[34],raw_bytes_from_sim[33] };
				//memcpy(&acc_z, &j, sizeof(acc_z));

				//acc_z = (raw_bytes_from_sim[33] << 24) | (raw_bytes_from_sim[34] << 16) | (raw_bytes_from_sim[35] <<8) | (raw_bytes_from_sim[36]) ;
				acc_z = conv_to_float(raw_bytes_from_sim[36], raw_bytes_from_sim[35],raw_bytes_from_sim[34],raw_bytes_from_sim[33]);


				//unsigned char k[] = { raw_bytes_from_sim[40], raw_bytes_from_sim[39],raw_bytes_from_sim[38],raw_bytes_from_sim[37] };
				//memcpy(&P, &k, sizeof(P));

				//P = (raw_bytes_from_sim[40] << 24) | (raw_bytes_from_sim[39] << 16) | (raw_bytes_from_sim[38] <<8) | (raw_bytes_from_sim[37]) ;
				P = conv_to_float(raw_bytes_from_sim[40], raw_bytes_from_sim[39],raw_bytes_from_sim[38],raw_bytes_from_sim[37]);

				//unsigned char l[] = { raw_bytes_from_sim[44], raw_bytes_from_sim[43],raw_bytes_from_sim[42],raw_bytes_from_sim[41] };
				//memcpy(&Q, &l, sizeof(Q));

				//Q = (raw_bytes_from_sim[41] << 24) | (raw_bytes_from_sim[42] << 16) | (raw_bytes_from_sim[43] <<8) | (raw_bytes_from_sim[44]) ;
				Q = conv_to_float(raw_bytes_from_sim[44], raw_bytes_from_sim[43],raw_bytes_from_sim[42],raw_bytes_from_sim[41]);

				//unsigned char m[] = { raw_bytes_from_sim[48], raw_bytes_from_sim[47],raw_bytes_from_sim[46],raw_bytes_from_sim[45] };
				//memcpy(&latitude, &m, sizeof(latitude));

				//latitude = (raw_bytes_from_sim[45] << 24) | (raw_bytes_from_sim[46] << 16) | (raw_bytes_from_sim[47] <<8) | (raw_bytes_from_sim[48]) ;
				latitude = conv_to_float(raw_bytes_from_sim[48], raw_bytes_from_sim[47],raw_bytes_from_sim[46],raw_bytes_from_sim[45]);

				//unsigned char n[] = { raw_bytes_from_sim[52], raw_bytes_from_sim[51],raw_bytes_from_sim[50],raw_bytes_from_sim[49] };
				//memcpy(&longitude, &n, sizeof(longitude));
				
				//longitude = (raw_bytes_from_sim[49] << 24) | (raw_bytes_from_sim[50] << 16) | (raw_bytes_from_sim[51] <<8) | (raw_bytes_from_sim[52]) ;
				longitude = conv_to_float(raw_bytes_from_sim[52], raw_bytes_from_sim[51],raw_bytes_from_sim[50],raw_bytes_from_sim[49]);


				
				// this will be paralellized
		        //Estimator( acc_x,  acc_y,  acc_z,  P, Q);

		        int slave_param =1;

				corethread_create(1, &Estimator1, (void*)slave_param); 
				corethread_create(2, &Estimator2, (void*)slave_param);

				corethread_join(1,(void*)slave_param);
				corethread_join(2,(void*)slave_param);



				//corethread_create(1, &Estimator, (void*)slave_param); 

				//corethread_join(1,(void*)slave_param);

				//Estimator();

				/*pitch_control_prop_gain = raw_bytes_from_sim[25];
				pitch_control_diff_gain = raw_bytes_from_sim[26];
				pitch_control_int_gain = raw_bytes_from_sim[27];

				pitch_control_prop_gain = pitch_control_prop_gain  / 10;
				pitch_control_diff_gain = pitch_control_diff_gain / 100;
				pitch_control_int_gain = pitch_control_int_gain / 10000;*/

				printf("Pitch: %f, Roll: %f, IAS: %f, accx: %f, accy: %f, accz: %f, P: %f Q: %f  lat: %f  lon: %f  thetac: %f  phic: %f\n",  
					pitch, roll, ias_sen, acc_x, acc_y, acc_z, P, Q, latitude, longitude, theta_cal, phi_cal);

				pitch = theta_cal;
				roll = phi_cal;

				//printf("%d %d %d %d \n",raw_bytes_from_sim[36], raw_bytes_from_sim[35], raw_bytes_from_sim[34], raw_bytes_from_sim[33]);

		
		}
		

		all_received = 0;



		if (mode == AUTO_TAKEOFF) {
			n_wheel = 0;
			throttle = 100;
			rudder = 0;
			elev = 05;
			ail = 0;
			if (ias_sen > 170) {
				mode = CLIMB;
			}
		}
		else if (mode == CLIMB) {
			if (alt_ft_msl <= 3000) {
				longitudinal_control(15);
				lateral_control(0);
			}
			else {
				mode = ALT_HOLD;
			}
		}
		else if (mode == ALT_HOLD) {
			HDG_control(150);
			cruise_control(4000);
		}

		


		elev = safe_byte(elev);
		ail = safe_byte(ail);
		rudder = safe_byte(rudder);
		throttle = safe_byte(throttle);

		/*out_elev_L = (elev & 0xFF);
		out_elev_H = ((elev >> 8) & 0xFF);

		out_ail_L = (ail & 0xFF);
		out_ail_H = ((ail >> 8) & 0xFF);

		out_rud_L = (rudder & 0xFF);
		out_rud_H = ((rudder >> 8) & 0xFF);

		out_thr_L = (throttle & 0xFF);
		out_thr_H = ((throttle >> 8) & 0xFF);*/

		BYTE send_sequence [] = { 255,254,253,252, elev, ail, rudder, throttle};

		//printf("Elev: %d,aileron: %d, rudder:  %d throttle: %d \n", elev, ail, rudder, throttle );


		for (int send_i = 0; send_i < 8; send_i ++){

             write(uart2_ptr,send_sequence[send_i]);
		}
		
		if(ias_sen < 20){
			mode = AUTO_TAKEOFF;
		}


	}//for
}
