/* Regularly spiking (RS) neuron according to Izhikevich
Nearest-neighbour, presynaptically centered, additive metaSTDP
copyright Lubica Benuskova
*/

// libraries used in the code 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

typedef struct _neuron *neuron;

neuron make_neuron_struct();
void init_neuron(neuron n);
neuron make_neuron();
void spontaneous_input(int time, neuron single_cell);
void cycle_cell(int time, neuron single_cell);
void stdp_train_cell(neuron single_cell);
void bcm_threshold(int time, neuron single_cell);
void test_sample(int time, neuron single_cell);
void run_network(int time);
void getTime(long *z);

#define SIMULATION_PERIOD  100000  // total # of time steps of the simulation in ms 

#define TIME_STEP               1    // time step in ms
#define SAMPLE_PERIOD       60000    // time step between sample points for plots
#define minutes             60000.0  // to get time in minutes in the format of real number

#define AP_peak                30.0   // peak of action potential for the given neuron  
#define P                       0.007 // probability P=0.00X corresponds to X Hz average frequency of homogenous Poisson series
#define MAX                 10000     // max number of presynaptic spikes between 2 postsynaptic spikes = 10000

#define TAU                 30000    // Integration period for theta_M in ms 
#define BOLTZMAN                0  	// should a Boltzman average be used
#define alpha                2000.0  // theta_M is scaled by this parameter 

// file variables for data used for plots
FILE *f_cell_voltage,*f_input_1,*f_input_2,*f_weight_1,*f_weight_2,*f_theta_M;

/**************************************************************************
		Variables
***************************************************************************/

/* weight initialization (Initial bias) */
float WEIGHT_1_INIT = 5;
float WEIGHT_2_INIT = 5;

/* STDP rule parameters */
float tau_ltp = 20;     // LTP tau in ms 
float tau_ltd = 100;    // LTD tau in ms 
float A0_ltp = 0.02;    // initial LTP peak at t=0
float A0_ltd = 0.01;    // initial LTD peak at t=0 

/* Izhikevich's neuron parameters values for RS cell*/
float coeff1 = 0.04;
float coeff2 = 5;
float coeff3 = 140;
float coeff4 = 1;

float a = 0.02;
float b = 0.2;
float c = -69.0;
float d = 2;

/* global variables */
int input_1 = 0;   // input = 0 if there is no presynaptic spike
int input_2 = 0;   // input = 1 if there is a presynaptic spike 
float rnd = 0.0;   // initial value of the random number for Poisson input
int t1_cell = 0;   // time of the previous postsynaptic spike
int t2_cell = 0;   // time of the following postsynaptic spike 
int index_1 = 0;   // cycle index for input 1
int index_2 = 0;   // cycle index for input 2
int t_input_1[MAX],t_input_2[MAX];  // arrays for times of presynaptic spikes 
float input_1_thetaM[MAX],input_2_thetaM[MAX]; // arrays for current values of theta_M 
int last_sample = 0;  // sampling of variables for plotting
float mean_out_frequency = 0.0; 

/**************************************************************************
		End of Variables and Parameters
***************************************************************************/

/***********************************************************************
		Neuron structure
***********************************************************************/

struct _neuron
{
/* variables describing the neuron  */
        float v;      // membrane voltage in mV 
        float u;      // membrane recovery variable
        int output;   // binary output
/* weights of neuron */
        float weight_1; // input 1 weight 
        float weight_2; // input 2 weight 
/* BCM averaging stuff */
        float theta_M;	// dynamic modification LTD-LTP threshold
        float mem;      // average of postsynaptic activity
        float pf; 	    // used to calculate mem
        float boltzfactor;	// used to calculate mem
}; /*end struct _neuron */

/************************************************************************
			Neuron Stuff
*************************************************************************/

/********* This function creates a new neuron and initializes it to 0 ******/
                   
neuron make_neuron_struct()
{
    neuron n;
 
    n = (neuron) malloc(sizeof(struct _neuron));

    n->v = c;
    n->u = b * c; 
    n->output = 0;	

    n->weight_1 = WEIGHT_1_INIT;
    n->weight_2 = WEIGHT_2_INIT; 
   
    n->theta_M = 1.0; 
    n->mem = 0.0;
    n->pf = 1.0;
    n->boltzfactor = 0.0; 

    return n;
} /*end neuron make_neuron_struct()*/

void init_neuron(neuron n)
{
    n->boltzfactor = exp( -((double) TIME_STEP / (double) TAU));
    if(BOLTZMAN == 0)
	n->pf = 1.0 / (1 - n->boltzfactor);
    else
	n->pf = 1.0;
} /* end void init_neuron(n)*/

/******** This function creates and initializes a new neuron **********/
neuron make_neuron()
{
    neuron n;

    n = make_neuron_struct();
    init_neuron(n);

    return n;
} /* end neuron make_neuron()*/

/************************************************************************
			End Neuron Stuff
*************************************************************************/

void main()
{ 
  long z;
  getTime(&z);
  srandom(z); // each time the code is run, it starts with a different seed for the random number generator, set according to real clock

  f_input_1 = fopen("input1_pulse_times.dat","w");
  f_input_2 = fopen("input2_pulse_times.dat","w"); 
  f_cell_voltage = fopen("cell_voltage.dat","w"); 
  f_weight_1 = fopen("weight_1.dat","w");
  f_weight_2 = fopen("weight_2.dat","w");
  f_theta_M = fopen("theta_M.dat","w");

  run_network(SIMULATION_PERIOD); // all computations happen here

  fclose(f_input_1);
  fclose(f_input_2);
  fclose(f_cell_voltage);
  fclose(f_weight_1);
  fclose(f_weight_2);
  fclose(f_theta_M);
      
} /* end void main()*/

/************************************************************************
			Network Stuff
*************************************************************************/

void run_network(int t)
{
    int time, i; 
        
    neuron single_cell;
    single_cell = make_neuron();
	
    // initialization of arrays 
    for(i = 0; i < MAX; ++i) 
    {t_input_1[i] = 0; t_input_2[i] = 0; input_1_thetaM[i] = 0.0; input_2_thetaM[i] = 0.0;}

    time = 0;
    while(time <= t)
    {
     spontaneous_input(time, single_cell);  // generation of inputs at time t
     cycle_cell(time, single_cell);         // update the cell according to the current inputs 
     
     if((time == t2_cell) && (t1_cell > 0)) stdp_train_cell(single_cell);  // train the cell according to the current inputs and output 

     bcm_threshold(time, single_cell);  // update theta_M for the next cycle
     test_sample(time, single_cell); // monitor variables 
	 
     time = time + TIME_STEP;      // increase time step by 1
    } /* while */ 

    mean_out_frequency = mean_out_frequency / ((float)SIMULATION_PERIOD / 1000.0); 
    printf("Average output frequency\t");
    printf("%f\n",mean_out_frequency); 
        
}  /* end  void run_network()*/ 

/************************************************************************
This function generates a spontaneous input
************************************************************************/
void spontaneous_input(int time, neuron n)
{
    float rnd; // random number generated from interval (0, 1)
     
    input_1 = 0; input_2 = 0; // initialization of inputs
  
// totally uncorrelated Poisson inputs 
    rnd = (float)(random())/(float)RAND_MAX;   
    if(P > rnd) 
       {
        input_1 = 1;  
        fprintf(f_input_1,"%f\t%d\n",(float)time,input_1); // data file to be filled with times of presynaptic spikes at input_1

        if(t1_cell > 0)
        { if((index_1 < MAX) && (time > t2_cell))  
        { t_input_1[index_1] = time; // array to be filled with times of presynaptic spikes at input_1
          input_1_thetaM[index_1] = n->theta_M; // array to be filled with values of theta_M when there is spike at input_1
          index_1 = index_1 + 1; 
        }
	    }
       }
                             
    rnd = (float)(random())/(float)RAND_MAX;   
    if(P > rnd) 
      {
        input_2 = 1;  
        fprintf(f_input_2,"%f\t%d\n",(float)time,input_2); // data file to be filled with times of presynaptic spikes at input_2

        if(t1_cell > 0)
        { if((index_2 < MAX) && (time > t2_cell))  
        { t_input_2[index_2] = time; // array to be filled with times of presynaptic spikes at input_2
          input_2_thetaM[index_2] = n->theta_M; // array to be filled with values of theta_M when there is spike at input_2
          index_2 = index_2 + 1; 
        }
	    }
       }              

/*
// correlated i.e. synchronous Poisson inputs
    rnd = (float)(random())/(float)RAND_MAX;   
    if(P > rnd) 
       {
        input_1 = 1; 
        input_2 = 1;
        fprintf(f_input_1,"%f\t%d\n",(float)time,input_1); // data file to be filled with times of presynaptic spikes at input_1
        fprintf(f_input_2,"%f\t%d\n",(float)time,input_2); // data file to be filled with times of presynaptic spikes at input_2

        if(t1_cell > 0)
        { if((index_1 < MAX) && (index_2 < MAX) && (time > t2_cell)) 
		{ t_input_1[index_1] = time;  // array to be filled with times of presynaptic spikes at input_1
	      t_input_2[index_2] = time;  // array to be filled with times of presynaptic spikes at input_2 
            
		  input_1_thetaM[index_1] = n->theta_M; // array to be filled with values of theta_M when there is spike at input_1
          input_2_thetaM[index_2] = n->theta_M; // array to be filled with values of theta_M when there is spike at input_2           
        
		  index_1 = index_1 + 1; 
          index_2 = index_2 + 1;
        }
	    }
      }
*/
              	
}  /* end void spontaneous_input(time) */


/**********************************************************************
This function updates the neuron output according to its inputs
**********************************************************************/
void cycle_cell(int time, neuron n)
{   
  n->output = 0;

   if(n->v >= AP_peak)  
     { if(t1_cell == 0) {t1_cell = time; t2_cell = time;} 
       else {t1_cell = t2_cell; t2_cell = time;} 
       n->v = c; n->u = n->u + d;
       mean_out_frequency = mean_out_frequency + 1;
     }
   else 
     {n->v = n->v + coeff1*(n->v * n->v) + coeff2*n->v + coeff3 - coeff4 * n->u;
      if(input_1 > 0) n->v = n->v + n->weight_1; 
      if(input_2 > 0) n->v = n->v + n->weight_2;   
      if(n->v >= AP_peak) {n->v = AP_peak; n->output = 1;}   // capping of output spike peak and setting the variable output to 1
      n->u = n->u + a*(b * n->v - n->u);  
     }
                                                                           
}  /* end void cycle_cell(time, n) */

/********************************************************************** 
                   STDP training
**********************************************************************/
void stdp_train_cell(neuron n)
{
  float ltp, ltd, A_ltp, A_ltd;
  int k;

  index_1 = 0; index_2 = 0; 

  while(index_1 < MAX)
    {if(t_input_1[index_1] == 0) index_1 = MAX;   // case when there are no presynaptic spikes at input_1
     else 
      {
       {
        if(input_1_thetaM[index_1]!=0) A_ltp = A0_ltp * (1.0/input_1_thetaM[index_1]); // update of LTP peak according to the current theta_M
        A_ltd = A0_ltd * (input_1_thetaM[index_1]); // update of LTD peak according to the current theta_M

        if(t2_cell > t_input_1[index_1])
        ltp = A_ltp * exp(-(t2_cell - t_input_1[index_1]) / tau_ltp); // calculation of actual increase in weight based on pre-post pairing
        if(t1_cell < t_input_1[index_1])
        ltd = A_ltd * exp((t1_cell - t_input_1[index_1]) / tau_ltd); // calculation of actual decrease in weight based on post-pre pairing

        n->weight_1 = n->weight_1 * (1 + ltp-ltd); // overall update of the weight 
		
		// the following lines may impose hard limits on weights if needed
        // if(n->weight_1 < 0.01) n->weight_1 = 0.01;
        // if(n->weight_1 > 0.09) n->weight_1 = 0.09;   
	
       }
        index_1 = index_1 + 1;
      }
    }

  while(index_2 < MAX)
    {if(t_input_2[index_2] == 0) index_2 = MAX; // case when there are no presynaptic spikes at input_2 
     else
      {     
       {
        if(input_2_thetaM[index_2]!=0) A_ltp = A0_ltp * (1.0/input_2_thetaM[index_2]);  
        A_ltd = A0_ltd * (input_2_thetaM[index_2]); 

        if(t2_cell > t_input_2[index_2])
        ltp = A_ltp * exp(-(t2_cell - t_input_2[index_2]) / tau_ltp); 
        if(t1_cell < t_input_2[index_2])
        ltd =  A_ltd * exp((t1_cell - t_input_2[index_2]) / tau_ltd);        

        n->weight_2 = n->weight_2 * (1 + ltp-ltd);	 	
        
		// the following lines may impose hard limits on weights if needed
		// if(n->weight_2 < 0.01) n->weight_2 = 0.01;
        // if(n->weight_2 > 0.09) n->weight_2 = 0.09;
		
        }
        index_2 = index_2 + 1; 
       } 
      }

  // resetting 
  for(k = 0; k < MAX; ++k) 
  {t_input_1[k] = 0; t_input_2[k] = 0; input_1_thetaM[k] = 0.0; input_2_thetaM[k] = 0.0;} 
  index_1 = index_2 = 0; 
  t1_cell = t2_cell; // important: time of the following postsynaptic spike become the time of the preceding postsynaptic spike

} /* end void stdp_train_cell(n) */


/********************************************************************** 
                   BCM theta_M calculation 
**********************************************************************/
void bcm_threshold(int time, neuron n)
{

  n->mem = (n->mem * (n->pf - 1.0)  + (n->output * n->output)) / n->pf;
  n->pf = n->pf * n->boltzfactor + 1.0;
  n->theta_M = alpha * n->mem;

} /* end void bcm_threshold(time, n) */

/************************************************************************
This function will test the cell every so often
************************************************************************/
void test_sample(int time, neuron n)
{ 
   
   if(n->output > 0) fprintf(f_cell_voltage,"%f\t%1d\n",(float)time/minutes,1);

  // saving data point for plots every minute, i.e. 60000 ms  
  if((time == 0) || (time - last_sample >= SAMPLE_PERIOD))
   {
    last_sample = time;    
    fprintf(f_weight_1,"%f\t%6.3f\n",(float)time/minutes,100*(n->weight_1-0.05)/0.05); // recaluculation as percents
    fprintf(f_weight_2,"%f\t%6.3f\n",(float)time/minutes,100*(n->weight_2-0.05)/0.05); // recaluculation as percents
    fprintf(f_theta_M,"%f\t%f\n",(float)time/minutes,n->theta_M); 	
   }

} /* end void test_sample(time, n) */

/***********************************************************************
This function returns the current real clock time in miniseconds for setting up the new seed of the random number generator
***********************************************************************/
void getTime(long *z) 
{
 struct timeval tp;
 struct timezone tzp;
  
 gettimeofday(&tp,&tzp);
 *z = tp.tv_usec;
}
    
