
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct particle //structure of a particle
{
	double mass; //particles mass
	double xlocation; //particles location on the x axis
	double ylocation; //particles location on the y axis
	double xvelocity; //particles velocity on the x axis
	double yvelocity; // particles velocity on the y axis
	double xacceleration; //particles acceleration on the x axis
	double yacceleration; //particles acceleration on the x axis
	int xwall; //checks whether particle has collided with left/right wall
	int ywall; //checks whether particles has collded with top/bottom wall
};

int getinfo(double *tini, double *tequi, double *tfin, double *dt, double *nsamp, double *contxl, double *contxr, double *contyb, double *contyt, double *epsilon, double *sigma, double *gravity, struct particle p[200])
{

//function that reads information from input files and saves them in the appropriate locations

int n = 0;

FILE *input; 
input = fopen("t7.txt", "r");

char dispensible[150]; //the dispensible string is used to temporarely store the titles of the input files, i.e. TIME, POTENTIAL

fscanf(input, "%s", dispensible); 
fscanf(input, "%lf", tini); //inital time
fscanf(input, "%lf", tequi); //equilibrium time
fscanf(input, "%lf", tfin); //final time
fscanf(input, "%lf", dt); //length of one timestep
fscanf(input, "%lf", nsamp); //number of samples

fscanf(input, "%s", dispensible);
fscanf(input, "%lf", contxl); //coordinates of left wall
fscanf(input, "%lf", contxr); //coordinates of right wall
fscanf(input, "%lf", contyb); //coordinates of bottom wall
fscanf(input, "%lf", contyt); //coordinates of top wall

fscanf(input, "%s", dispensible);
fscanf(input, "%lf", epsilon); 
fscanf(input, "%lf", sigma);
fscanf(input, "%lf", gravity);

fscanf(input, "%s", dispensible);

while (!feof(input)) //stores information about each particle in the struct defined earlier, loop keeps going untill it reached the end of the file
{
	fscanf(input, "%lf", &p[n].mass);
	fscanf(input, "%lf", &p[n].xlocation);
	fscanf(input, "%lf", &p[n].ylocation);
	fscanf(input, "%lf", &p[n].xvelocity);
	fscanf(input, "%lf", &p[n].yvelocity);
	n++;
}

fclose(input);


return(n-1);


}


int checkerror(int particlenum, double *tini, double *tequi, double *tfin, double *dt, double *nsamp, double *contxl, double *contxr, double *contyb, double *contyt, double *epsilon, double *sigma, struct particle p[200])
{
	int n=0;

	//function that checks whether there are any errors in the file

	//this part checks errors not regarding the particles. the errors are as follows: final time smaller than initial or equilbrium times, timestep duration shorter than zero, negative number of samples, negative box size, negative sigma (has units of length)
	if(((*tfin)<(*tini))||((*tequi)<(*tini))||((*dt)<=0)||((*contxr)<(*contxl))||((*contyt)<(*contyb))||((*sigma)<0)||((*nsamp)<=0))
		{
			n=1;
		}

	//this part checks errors regarding the particles. the errors are as follows: negative mass, particles initial location outside of box
	for (int i = 0; i < particlenum; i++)
	{
		if(p[i].mass<0||p[i].xlocation<(*contxl)||p[i].xlocation>(*contxr)||p[i].ylocation<(*contyb)||p[i].ylocation>(*contyt))
		{
			n++;
		}
	}

	if(n==0)
	{
		return 0;
	}
	else
		return 1;
}

double checkkinetic (struct particle p[200], int particlenum, double dt)
{

	//function that checks the average kinetic energy of all particles in the container every timestep and returns it to a double in the main function

	double k=0;
	for (int i = 0; i < particlenum; i++)
	{	
		k=k+p[i].mass*(pow(p[i].xvelocity+p[i].xacceleration*(dt/2),2)+pow(p[i].yvelocity+p[i].yacceleration*(dt/2),2))/2; 
	}
	return (k/particlenum);
}

double checkpotential (struct particle p[200], int particlenum, double *contyb, double *gravity, double *sigma, double *epsilon)
{

	//function that checks the average potential energy of all particles in the container every timestep and returns it to a double in the main function

	double u=0, h=0, d;

	//potential from the lennard-jones potential
	for (int i = 0; i < particlenum; i++)
	{
		for(int j = i; j < particlenum ; j++)
		{
			if((i!=j)&&(d<3.54938*(*sigma)))
			{
				d=sqrt(pow(p[i].xlocation-p[j].xlocation,2)+pow(p[i].ylocation-p[j].ylocation,2));
				u=u+((*epsilon)*((pow((*sigma),12)/pow(d,12))-2*(pow((*sigma),6)/pow(d,6))+pow(10,-3)));
			}
		}

		//potential from gravity
		h=h+p[i].mass*(*gravity)*(p[i].ylocation-(*contyb));
	}
	return ((h+u)/particlenum);
}

double checkpress (struct particle p[200], int particlenum, double *dt, double length)
{
	//function that checks the pressure applied to the conatiners walls every timestep

	double t=0;

	for (int i = 0; i < particlenum; i++)
	{
		if(p[i].xwall==1) //gets input on whether the particle has collided with a wall
		{
			t = t + (((p[i].xvelocity)*(p[i].mass))*2)/((*dt)*length);
		}
		if(p[i].ywall==1) //gets input on whether the particle has collided with a wall
		{
			t = t + (((p[i].yvelocity)*(p[i].mass))*2)/((*dt)*length);
		}
	}

	return t;
}


void xacceleration (struct particle p[200], int particlenum, double *sigma, double *epsilon)
{
	
	//function that calculates the acceleration of each particle in the x axis according to the deriviative of the lennard jones potential

	double d ;
	for (int i = 0; i < particlenum; i++)
	{		
		for (int j = 0; j < particlenum; j++)
		{
			d=sqrt(pow((p[i].xlocation-p[j].xlocation),2)+pow((p[i].ylocation-p[j].ylocation),2));
			if((i!=j)&&(d<3.54938*(*sigma)))
			{
				p[i].xacceleration+=((12*(*epsilon)*((pow(*sigma,12))/(pow(d,14))-((pow(*sigma,6)/(pow(d,8))))))*(p[i].xlocation-p[j].xlocation))/p[i].mass;

			}
		}
	}
}	

void yacceleration (struct particle p[200], int particlenum, double *sigma, double *epsilon, double *gravity)
{

	//function that calculates the acceleration of each particle in the y axis according to the deriviative of the lennard jones potential

	double d ;
	for (int i = 0; i < particlenum; i++)
	{
		for (int j = 0; j < particlenum; j++)
		{
			d=sqrt(pow(p[i].xlocation-p[j].xlocation,2)+pow(p[i].ylocation-p[j].ylocation,2));
			if((i!=j)&&(d<3.54938*(*sigma)))
			{
				p[i].yacceleration+=((12*(*epsilon)*((pow((*sigma),12))/(pow(d,14))-((pow((*sigma),6)/(pow(d,8))))))*(p[i].ylocation-p[j].ylocation))/p[i].mass;	
			}
		}
	}
}
		
void xparlocation(struct particle p[200], int particlenum, double *dt)
{

	//function that calculates the location of each particle in the x axis according to its former location and velocity

	for (int i = 0; i < particlenum; i++)
	{
		p[i].xlocation=p[i].xlocation+p[i].xvelocity*(*dt);
	}
}

void yparlocation(struct particle p[200], int particlenum, double *dt)	
{

	//function that calculates the location of each particle in the y axis according to its former location and velocity

	for (int i = 0; i < particlenum; i++)
	{
		p[i].ylocation=p[i].ylocation+p[i].yvelocity*(*dt);
	}
}	

void inixvelocity(struct particle p[200], int particlenum, double *dt)
{

	//function that calculates the velocity of each particle in the x axis according to its former velocity and acceleration

	for (int i = 0; i < particlenum; i++)
	{
		p[i].xvelocity=p[i].xvelocity+p[i].xacceleration*((*dt)/2);
	}
}

void iniyvelocity(struct particle p[200], int particlenum, double *dt)
{

	//function that calculates the velocity of each particle in the y axis according to its former velocity and acceleration

	for (int i = 0; i < particlenum; i++)
	{
		p[i].yvelocity=p[i].yvelocity+p[i].yacceleration*((*dt)/2);
	}
}

void xvelocity(struct particle p[200], int particlenum, double *dt)
{

	//function that calculates the velocity of each particle in the x axis according to its former velocity and acceleration

	for (int i = 0; i < particlenum; i++)
	{
		p[i].xvelocity=p[i].xvelocity+p[i].xacceleration*(*dt);
	}
}

void yvelocity(struct particle p[200], int particlenum, double *dt)
{

	//function that calculates the velocity of each particle in the y axis according to its former velocity and acceleration

	for (int i = 0; i < particlenum; i++)
	{
		p[i].yvelocity=p[i].yvelocity+p[i].yacceleration*(*dt);
	}
}

void wallcollision (struct particle p[200], int particlenum, double *contxl, double *contxr, double *contyb, double *contyt)
{

	//function that checks whether a particle has collided with a wall
	//if function detects a collision, it changes its velocity to the opposite direction, updates its location to other side of wall (with same distance from wall), and updates an integer that indicates that a collision has occured

	for (int i = 0; i < particlenum; i++)
	{
		if(p[i].xlocation<=(*contxl))
		{
			p[i].xvelocity=(-1)*p[i].xvelocity;
			p[i].xlocation=2*(*contxl)-p[i].xlocation;
			p[i].xwall=1;
		}
		if(p[i].xlocation>=(*contxr))
		{
			p[i].xvelocity=(-1)*p[i].xvelocity;
			p[i].xlocation=2*(*contxr)-p[i].xlocation;
			p[i].xwall=1;
		}
		if (p[i].ylocation<=(*contyb))
		{
			p[i].yvelocity=(-1)*p[i].yvelocity;
			p[i].ylocation=2*(*contyb)-p[i].ylocation;
			p[i].ywall=1;
		}
		if (p[i].ylocation>=(*contyt))
		{
			p[i].yvelocity=(-1)*p[i].yvelocity;
			p[i].ylocation=2*(*contyt)-p[i].ylocation;
			p[i].ywall=1;
		}
	}
}

void simulateout(double i, int sample, int particlenum, struct particle p[200], double *nsamp, double *contxl, double *contxr, double *contyb, double *contyt)
{

	//function that creates and output file

	FILE *output; 
	output = fopen("C:\\Users\\Julian\\Desktop\\simulate.dat", "a");

	if (sample==0)
	{
		fprintf(output,"%lf", (*nsamp));
		fprintf(output,"	%lf", (*contxl));
		fprintf(output,"	%lf", (*contxr)); 
		fprintf(output,"	%lf", (*contyb));
		fprintf(output,"	%lf", (*contyt));
		fprintf(output,"	%d \r\n \r\n", particlenum);
	}

	fprintf(output,"%d", sample);
	fprintf(output,"	%lf \r\n", i);

	for (int n = 0; n < particlenum; n++)
	{
		fprintf(output,"%lf", p[n].mass);
		fprintf(output,"	%lf", p[n].xlocation); 
		fprintf(output,"	%lf", p[n].ylocation);
		fprintf(output,"	%lf", p[n].xvelocity);
		fprintf(output,"	%lf \r\n", p[n].yvelocity);
	}

	fprintf(output,"\r\n");

}


void testout(double i, int sample, int particlenum, struct particle p[200], double *nsamp, double *contxl, double *contxr, double *contyb, double *contyt)
{

	//function that creates and output file

	FILE *output; 
	output = fopen("C:\\Users\\Julian\\Desktop\\test.dat", "a");

	if(particlenum<5)
	{
		if(sample==0)
		{
		fprintf(output, "4 particles or less, file irrelevant. see simulate file.");
		}
	}	

	else
	{

		if (sample==0)
	{
		fprintf(output,"%lf", (*nsamp));
		fprintf(output,"	%lf", (*contxl));
		fprintf(output,"	%lf", (*contxr)); 
		fprintf(output,"	%lf", (*contyb));
		fprintf(output,"	%lf", (*contyt));
		fprintf(output,"	%d \r\n \r\n", particlenum);
	}

	fprintf(output,"%d", sample);
	fprintf(output,"	%lf \r\n", i);

	for (int n = 0; n < 2; n++)
	{
		fprintf(output,"%lf", p[n].mass);
		fprintf(output,"	%lf", p[n].xlocation); 
		fprintf(output,"	%lf", p[n].ylocation);
		fprintf(output,"	%lf", p[n].xvelocity);
		fprintf(output,"	%lf \r\n", p[n].yvelocity);
	}

	for (int n = (particlenum-2); n < particlenum; n++)
	{
		fprintf(output,"%lf", p[n].mass);
		fprintf(output,"	%lf", p[n].xlocation); 
		fprintf(output,"	%lf", p[n].ylocation);
		fprintf(output,"	%lf", p[n].xvelocity);
		fprintf(output,"	%lf \r\n", p[n].yvelocity);
	}

	fprintf(output,"\r\n");

	}

}


void resultout(double i, int sample, double avgtemp, double avgpres, double ratio1, double ratio2)
{

	//function that creates and output file

	FILE *output; 
	output = fopen("C:\\Users\\Julian\\Desktop\\results.dat", "a");

	fprintf(output,"%d", sample);
	fprintf(output,"	%lf", i);
	fprintf(output,"	%lf", avgtemp);
	fprintf(output,"	%lf", avgpres);
	fprintf(output,"	%lf", ratio1);
	fprintf(output,"	%lf\n", ratio2);
}


void energyout(double i, int sample, double ekinetic, double epotential, double etotal)
{

	//function that creates and output file

	FILE *output; 
	output = fopen("C:\\Users\\Julian\\Desktop\\energies.dat", "a");

	fprintf(output,"%d", sample);
	fprintf(output,"	%lf", i); 
	fprintf(output,"	%lf", ekinetic);
	fprintf(output,"	%lf", epotential);
	fprintf(output,"	%lf \r\n", etotal);

}


void momentaout(double i, int sample, int particlenum, struct particle p[200])
{

	//function that creates and output file

	FILE *output; 
	output = fopen("C:\\Users\\Julian\\Desktop\\momenta.dat", "a");

	double px=0, py=0;

	for (int n = 0; n < particlenum; n++)
	{
		px=px+((p[n].mass)*(p[n].xvelocity));
	}

	for (int n = 0; n < particlenum; n++)
	{
		py=py+((p[n].mass)*(p[n].yvelocity));
	}

	fprintf(output,"%d", sample);
	fprintf(output,"	%lf", i);
	fprintf(output,"	%lf", px);  
	fprintf(output,"	%lf \r\n", py);
}


int main()
{	

	//this part of the main is where the initial paramaters are defined and info from an input file gets stored in them
	double tini, tequi, tfin, dt, nsamp, contxl, contxr, contyb, contyt, epsilon, sigma, gravity;
	struct particle p[200];  
	int particlenum = getinfo(&tini, &tequi, &tfin, &dt, &nsamp, &contxl, &contxr, &contyb, &contyt, &epsilon, &sigma, &gravity, p);	
	int error = checkerror(particlenum, &tini, &tequi, &tfin, &dt, &nsamp, &contxl, &contxr, &contyb, &contyt, &epsilon, &sigma, p);

	//this part of the main is where the total length of the walls, the area of the container, the ratios, the parameters and the avarage temperature and pressure are defined	
	double length = (contxr-contxl)*2 + (contyt-contyb)*2, area = (contxr-contxl)*(contyt-contyb), ratio1, ratio2, a1, b1, avgtemp=0, avgpres=0;

	//this part calls for the inital calculation of the pressure, kinetic , potential and total energy of the system
	double press = checkpress (p, particlenum, &dt, length);
	double ekinetic = checkkinetic(p, particlenum, dt);
	double epotential = checkpotential(p, particlenum, &contyb, &gravity, &sigma, &epsilon);
	double etotal = ekinetic + epotential;

	//definitions of parameters of later use	
	//k-number of timesteps in this run, timestep-amount of time between each timestep, i-current time, sample-current sample
	//n,c,t are counters for different purposes within the loop
	int n=0, k=((tfin-tini)/dt), timestep=k/(nsamp-1), t=0, sample=0, c=0;
	double i=tini;

	
	//calls for printing of the first line of output files with the initial conditions
	simulateout(i, sample, particlenum, p, &nsamp, &contxl, &contxr, &contyb, &contyt);
	testout(i, sample, particlenum, p, &nsamp, &contxl, &contxr, &contyb, &contyt);	
	energyout(i, sample, ekinetic, epotential, etotal);
	momentaout(i, sample, particlenum, p);	

	//this calls for the initial calculation of the acceleration according to the initial conditions
	xacceleration(p, particlenum, &sigma, &epsilon);
	yacceleration(p, particlenum, &sigma, &epsilon, &gravity);

	//this function updates the velocity to that of dt/2 according to the initial acceleration and velocity, as done in the verlet leap from method
	inixvelocity(p, particlenum, &dt);
	iniyvelocity(p, particlenum, &dt);
	
	//loop that starts the movment of the system from the initial time + dt to the final time on intervals of dt, on the condition that there are no errors in the input
	for (i = (tini+dt); (i < tfin)&&(error==0) ; i+=dt)
	{		

			//calls for the calculation of the new location according to the velocity of dt/2 ago
			xparlocation(p, particlenum, &dt);
			yparlocation(p, particlenum, &dt);

			//these loops zero the acceleration each dt, so it is recalculated from scratch each timestep
			for (int j = 0; j < particlenum; j++)
			{
				p[j].xacceleration=0;
			}
			for (int j = 0; j < particlenum; j++)
			{
				p[j].yacceleration=(-1)*(gravity);
			}

			//this part of the calculation calls for the calculation of the acceleration and velocity of each particle each dt, as well as checks whether a wall collision has occured and wether the equilbrium time has arrived
			xacceleration(p, particlenum, &sigma, &epsilon);
			yacceleration(p, particlenum, &sigma, &epsilon, &gravity);

			xvelocity(p, particlenum, &dt);
			yvelocity(p, particlenum, &dt);
			
			wallcollision(p, particlenum, &contxl, &contxr, &contyb, &contyt);
			press = checkpress (p, particlenum, &dt, length);

			//this if condition checks whether the time has reached time equilbrium, and if so, begins counting the average temperature and pressure of the system 
		    if(i>=tequi)
		    {
				c++;
				avgtemp = (avgtemp*(c-1) + ekinetic)/c;
				avgpres = (avgpres*(c-1) + press)/c;
				a1=1.4027*epsilon*pow(sigma,2)/(cbrt(epsilon/avgtemp)); 
				b1=1.7731*pow(sigma,2)*sqrt(cbrt(epsilon/avgtemp));
				ratio1 = (avgpres*area)/(particlenum*avgtemp);
				ratio2 = (avgpres + (a1*pow(particlenum,2))/(pow(area,2)))*(area-(particlenum*b1))/(particlenum*avgtemp);
		    }		
		t++;

		//calls for the calculation of the kinetic and potential energy each dt/2 to keep it updated
		ekinetic = checkkinetic(p, particlenum, dt);
		epotential = checkpotential(p, particlenum, &contyb, &gravity, &sigma, &epsilon);
		etotal = ekinetic + epotential;
		
		//this if condition checks whether it is time to print out a sample, according to the run time and number of samples requested
		if(t==timestep)
		{
			sample++;

			simulateout(i, sample, particlenum, p, &nsamp, &contxl, &contxr, &contyb, &contyt);
			testout(i, sample, particlenum, p, &nsamp, &contxl, &contxr, &contyb, &contyt);	
			energyout(i, sample, ekinetic, epotential, etotal);
			momentaout(i, sample, particlenum, p);
			
			//this if condition checks whether equilbrium time has arrived, and if so, starts printing out to an output files information regarding pressure, temperature and ratios
			if(i>=tequi)
			{		
			resultout(i, sample, avgtemp, avgpres, ratio1, ratio2);
			}
			
			t=0;
		}

		//this loop makes all integers indicating wall collision return to zero
		for (int j = 0; j < particlenum; j++)
		{
			p[j].xwall=0;
			p[j].ywall=0;
		}

	}

	//this if condition checks whether there was an error in the input files, and if so, returns and input to the user that indicates so
	if(error==1)
	{
		printf("error in input file");
	}
}


