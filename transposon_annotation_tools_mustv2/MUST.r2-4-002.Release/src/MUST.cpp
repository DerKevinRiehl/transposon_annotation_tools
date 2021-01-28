#include <fstream>
#include <iostream>
#include <cmath>
#include <string.h>
#include <stdlib.h>

#define MaxLineLength 10000

using namespace std;

//#define ERROR_number 1  ///use the error rate to calculate
//const int MIN_TIR_length= 9;
//const int MAX_TIR_length= 30;
//const int MIN_DR_length= 2;
//const int MAX_DR_length= 30;   
//const int MIN_IR_length= 60;
//const int MAX_IR_length= 600;

//const int FIXED_FLANKING_length=50;

int main(int argc, char *argv[])
{
	char filename[MaxLineLength];
	char gSeqID[MaxLineLength];
	char fout[MaxLineLength];
	char ftir[MaxLineLength];
	int i,j,k,m;
	int count_number_1=0;
	int count_number_2=0;
	int real_DR_length=0;
	int real_TIR_length=0;
	int mite_ID=0;
	double DR_correct_rate=0;
	double TIR_correct_rate=0;
	double DR_AT_rate=0;
	double TIR_AT_rate=0;
	double IR_AT_rate=0;
	double DRF_AT_rate=0;

	int stringlength;
	char * string_record;   

	//Initialization
	//if( (argc!=11) && (argc!=4) )
	if( (argc!=13) )
	{
		printf("Error in syntax!\n");
		printf("    ./MUST <input_seq.raw_data> <SeqID> <prediction.dat> <tir_seq.dat> [7 options]\n");
        printf("Option list:\n");
        printf("       <MIN_TIR_length>\n");
        printf("       <MAX_TIR_length>\n");
        printf("       <MIN_DR_length>\n");
        printf("       <MAX_DR_length>\n");
/*
        printf("       <MIN_IR_length>\n");
        printf("       <MAX_IR_length>\n");
*/
        printf("       <MIN_MITE_length>\n");
        printf("       <MAX_MITE_length>\n");
        printf("       <FIXED_FLANKING_length>\n");
		printf("       <MUTATION_RATE>\n");
		exit(-1);
	}
	
	strncpy(filename, argv[1], MaxLineLength);
	strncpy(gSeqID, argv[2], MaxLineLength);
	strncpy(fout, argv[3], MaxLineLength);
	strncpy(ftir, argv[4], MaxLineLength);

        int MIN_TIR_length=atoi(argv[5]);
        int MAX_TIR_length=atoi(argv[6]);
        int MIN_DR_length=atoi(argv[7]);
        int MAX_DR_length=atoi(argv[8]);
        int MIN_IR_length=atoi(argv[9])-MIN_TIR_length;
        int MAX_IR_length=atoi(argv[10])-MIN_TIR_length;
        int FIXED_FLANKING_length=atoi(argv[11]);
		double FIXED_ERROR_rate = atof(argv[12]);

	ifstream infile(filename,ios::in);
	if(!infile)
	{
		printf("can't open!",filename);
		return 0;
	}
	else
	{   char cc;	
		i=0;
		while(infile.get(cc))
		{
		   if(cc!='\n')
	    	i++;
		}
		stringlength=i;         //
	}
	cout<<endl;
	cout<<stringlength<<endl;

	infile.close();
	string_record= new char[3*FIXED_FLANKING_length+stringlength];	

	ifstream infile1(filename,ios::in);	

	for(i=0;i<FIXED_FLANKING_length;i++)  ///to first add FIXED_FLANKING_length number '-' in the haed of string_record
		string_record[i]='-';
		

	char ccc;
	i=FIXED_FLANKING_length;
	while(infile1.get(ccc))
	{
	    if(ccc!='\n')
		{
			string_record[i]=ccc;              //
			i++;
		}
	}

	infile1.close();
	//ofstream outfile("optimized_to_follow_text_temp.txt",ios::out);
	ofstream outfile(fout,ios::out);
	ofstream tofollow(ftir,ios::out);
        for(i=FIXED_FLANKING_length+stringlength;i<3*FIXED_FLANKING_length+stringlength;i++)  ///to  add 2*FIXED_FLANKING_length number '-' in the end of string_record
                string_record[i]='-';
        int left_position_controller;
        int right_position_controller;
	for(i=FIXED_FLANKING_length+MIN_DR_length;i<stringlength-MIN_IR_length-MIN_TIR_length-MIN_DR_length+FIXED_FLANKING_length;i++)  //make sure to get the 50 flanking 
	{
		if(i<FIXED_FLANKING_length+MAX_DR_length)
		{ left_position_controller=i-FIXED_FLANKING_length; }
		else { left_position_controller=MAX_DR_length; }
                if(i<FIXED_FLANKING_length+stringlength-MIN_DR_length-MIN_TIR_length-MAX_IR_length)
                //{ right_position_controller=MAX_IR_length; }
                { right_position_controller=MAX_IR_length-MIN_IR_length; }
                else{ right_position_controller=FIXED_FLANKING_length+stringlength-MIN_DR_length-MIN_TIR_length-MIN_IR_length-i;}

		real_TIR_length=MIN_TIR_length;
		DR_correct_rate=0;
		TIR_correct_rate=0;

		double FIXED_ERROR_rate_DR = FIXED_ERROR_rate;
		FIXED_ERROR_rate_DR = 1.0;
		for(k=i+MIN_IR_length;k<i+MIN_IR_length+right_position_controller;k++)
		{
	                real_DR_length=0;
			count_number_1=0;
			for(m=0;m<MIN_TIR_length;m++)
			{
				if((string_record[i+m]=='A')&&(string_record[k+MIN_TIR_length-m-1]=='T'))
					count_number_1++;
                else if((string_record[i+m]=='T')&&(string_record[k+MIN_TIR_length-m-1]=='A'))
                    count_number_1++;
                else if((string_record[i+m]=='C')&&(string_record[k+MIN_TIR_length-m-1]=='G'))
                    count_number_1++;
                else if((string_record[i+m]=='G')&&(string_record[k+MIN_TIR_length-m-1]=='C'))
                    count_number_1++;
			}

			// The boundary nucleotide must match!
			m=0;
				if((string_record[i+m]=='A')&&(string_record[k+MIN_TIR_length-m-1]!='T'))
					count_number_1 = 0;
                else if((string_record[i+m]=='T')&&(string_record[k+MIN_TIR_length-m-1]!='A'))
                    count_number_1 = 0;
                else if((string_record[i+m]=='C')&&(string_record[k+MIN_TIR_length-m-1]!='G'))
                    count_number_1 = 0;
                else if((string_record[i+m]=='G')&&(string_record[k+MIN_TIR_length-m-1]!='C'))
                    count_number_1 = 0;

			m=1;
				if((string_record[i+m]=='A')&&(string_record[k+MIN_TIR_length-m-1]!='T'))
					count_number_1 = 0;
                else if((string_record[i+m]=='T')&&(string_record[k+MIN_TIR_length-m-1]!='A'))
                    count_number_1 = 0;
                else if((string_record[i+m]=='C')&&(string_record[k+MIN_TIR_length-m-1]!='G'))
                    count_number_1 = 0;
                else if((string_record[i+m]=='G')&&(string_record[k+MIN_TIR_length-m-1]!='C'))
                    count_number_1 = 0;

			//if(count_number_1>=MIN_TIR_length-floor(MIN_TIR_length*0.1))
   	        if(count_number_1>=MIN_TIR_length*FIXED_ERROR_rate)
			{
				real_DR_length=0;
				for(j=MIN_DR_length;j<left_position_controller+1;j++)
				{
	                count_number_2=0;
    	            for(m=0;m<j;m++)
        	        {
            	    	if((string_record[i-j+m]==string_record[k+m+MIN_TIR_length])&&(string_record[i-j+m]!='N'))
                	    	count_number_2++;
                	}
					//if(count_number_2>=j-floor(j*0.1))

					//if(count_number_2>=j*FIXED_ERROR_rate)
					if(count_number_2>=j*FIXED_ERROR_rate_DR)
					{
						real_DR_length=j;
						DR_correct_rate=(double)count_number_2/(double)(j);
                        //TIR_correct_rate=(double)count_number_1/(double)(MIN_TIR_length);
						TIR_correct_rate=count_number_1;
				  	}
				}
			}

			if(real_DR_length>=MIN_DR_length)
			{
				count_number_1=0;
			
				for(j=MIN_TIR_length;j<MAX_TIR_length;j++)
				{
					if((string_record[i+j]=='A')&&(string_record[k-(j-MIN_TIR_length)-1]=='T'))
						continue;
	                                else if((string_record[i+j]=='C')&&(string_record[k-(j-MIN_TIR_length)-1]=='G'))
        	                                continue;
                	                else if((string_record[i+j]=='G')&&(string_record[k-(j-MIN_TIR_length)-1]=='C'))
                        	                continue;
                                	else if((string_record[i+j]=='T')&&(string_record[k-(j-MIN_TIR_length)-1]=='A'))
                                        	continue;
					else {count_number_1++;}
					if(count_number_1>=1)
					{
						real_TIR_length=j;
						break;
					}	
				}

				TIR_correct_rate=(double)(TIR_correct_rate+real_TIR_length-MIN_TIR_length)/(double)(real_TIR_length);
				
				count_number_1=0;
				for(j=i-real_DR_length;j<i;j++)
				{
					if((string_record[j]=='A')||(string_record[j]=='T'))
						count_number_1++;
				}
				DR_AT_rate=(double)count_number_1/(double)real_DR_length;
        	                count_number_1=0;
				for(j=i;j<i+real_TIR_length;j++)
				{
					if((string_record[j]=='A')||(string_record[j]=='T'))
                                        	count_number_1++;
				}
	                        TIR_AT_rate=(double)count_number_1/(double)real_TIR_length;
        	                count_number_1=0;
                	        for(j=i;j<k+MIN_TIR_length;j++)
                        	{
                                	if((string_record[j]=='A')||(string_record[j]=='T'))
                                        	count_number_1++;
                        	}
	                        IR_AT_rate=(double)count_number_1/(double)(k+MIN_TIR_length-i);
				count_number_1=0;
                       		for(j=i-real_DR_length-FIXED_FLANKING_length;j<i;j++)
                       		{
                               		if((string_record[j]=='A')||(string_record[j]=='T'))
                                       		count_number_1++;
                       		}
				DRF_AT_rate=(double)count_number_1/(double)(real_DR_length+FIXED_FLANKING_length);

				outfile<<"GenomeFile\t"<<filename<<"\nID\t"<<mite_ID << endl;
				outfile<<"GID\t" << gSeqID << endl;
				mite_ID++;
				outfile<<"POSITION\t"<<i-FIXED_FLANKING_length+1<<"-"<<k+MIN_TIR_length-FIXED_FLANKING_length << endl;  // +1 is to be the same in the NCBI
				outfile<<"LENGTH\t"<<abs((k+MIN_TIR_length-FIXED_FLANKING_length)-(i-FIXED_FLANKING_length+1))+1 << endl;
				outfile<<"DR\t"<<real_DR_length << endl;
				outfile<<"TIR\t"<<real_TIR_length << endl;
				outfile<<"DRIDENT\t"<<DR_correct_rate<<"\nTIRIDENT\t"<<TIR_correct_rate << endl;
				outfile<<"DRLEFT\t";
				for(m=0;m<real_DR_length;m++)
					outfile<<string_record[i-real_DR_length+m];
                	        outfile<<"\nTIRLEFT\t";
                        	for(m=0;m<real_TIR_length;m++)
                                {
					outfile<<string_record[i+m];
                                        tofollow<<string_record[i+m];
				}

				tofollow<<endl;
				outfile << endl;
				outfile<<"TIRRIGHT\t";
				for(m=0;m<real_TIR_length;m++)
					outfile<<string_record[k+m-(real_TIR_length-MIN_TIR_length)];
                        	outfile<<"\nDRRIGHT\t";
	                        for(m=0;m<real_DR_length;m++)
        	                        outfile<<string_record[k+m+MIN_TIR_length];
				outfile<<"\nDRAT\t"<<DR_AT_rate<<"\nTIRAT\t"<<TIR_AT_rate;
				
				float MITEAT = 0;
                for(j=i;j<k+MIN_TIR_length;j++)
       	        {
                    if( (string_record[j]=='A') || (string_record[j]=='T') )
                    {
                        MITEAT++;
                    }
                }
                MITEAT = MITEAT/(abs(k+MIN_TIR_length-i)+1);
				
				//outfile<<"\nIRAT\t"<<IR_AT_rate<<"\nFDRAT\t"<<DRF_AT_rate;
				outfile<<"\nMITEAT\t"<<MITEAT<<"\nFDRAT\t"<<DRF_AT_rate;
				outfile<<"\nPRESITE\t";
        	                for(j=i-real_DR_length-FIXED_FLANKING_length;j<i;j++)
                	        {
                        	      outfile<<string_record[j];
                        	}
	                        for(j=k+MIN_TIR_length+real_DR_length;j<k+MIN_TIR_length+real_DR_length+FIXED_FLANKING_length;j++)
        	                {
                	              outfile<<string_record[j];
                        	}
	                        outfile<<"\nPOSTSITE\t";
        	                for(j=i-real_DR_length-FIXED_FLANKING_length;j<i;j++)
                	        {
                        	      outfile<<string_record[j];
                        	}
	                        for(j=k+MIN_TIR_length;j<k+MIN_TIR_length+real_DR_length+FIXED_FLANKING_length;j++)
        	                {
                	              outfile<<string_record[j];
                        	}

				outfile<<"\nMITE\t";
        	                for(j=i;j<k+MIN_TIR_length;j++)
                	        {
                        	        outfile<<string_record[j];
                        	}
				outfile<<"\n//\n\n";
			}	
		}
	}
	outfile.close();
	tofollow.close();
	delete [] string_record;

return 0;
}
	
	
	
	
