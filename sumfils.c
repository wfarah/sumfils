#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mysigproc/mysigproc_utils.h"

#define NSEC 40.
#define MAXFILES 64

#ifndef MIN
#define MIN(a,b) ((a)<(b) ?(a):(b))
#endif

void usage()
{
  fprintf(stdout,
      "Usage:\n"
      "sumfils nfiles -o outfile\n"
      );
}

int sanity_check_files(fil_t ** infilbanks, int ninfiles)
{
  float fch1 = infilbanks[0]->header->fch1;
  float foff = infilbanks[0]->header->foff;
  float tsamp = infilbanks[0]->header->tsamp;
  float tstart = infilbanks[0]->header->tstart;

  for (int i=1; i<ninfiles; i++)
  {
    if (fch1 != infilbanks[i]->header->fch1)
    {
      fprintf(stderr, "ERROR: fch1 is not the same in all filbanks "
          "(%f, %f) \n", fch1, infilbanks[i]->header->fch1);
      return EXIT_FAILURE;
    }
    if (foff != infilbanks[i]->header->foff)
    {
      fprintf(stderr, "ERROR: foff is not the same in all filbanks "
          "(%f, %f) \n", foff, infilbanks[i]->header->foff);
      return EXIT_FAILURE;
    }
    if (abs(tsamp - infilbanks[i]->header->tsamp) > 1e-9)
    {
      fprintf(stderr, "ERROR: tsamp is not the same in all filbanks "
          "(%f, %f) \n", tsamp, infilbanks[i]->header->tsamp);
      return EXIT_FAILURE;
    }
    if (abs(tstart - infilbanks[i]->header->tstart) > 1e-9)
    {
      fprintf(stderr, "ERROR: tstart is not the same in all filbanks "
          "(%f, %f) \n", tstart, infilbanks[i]->header->tstart);
      return EXIT_FAILURE;
    }
  }
}


// 8bit
static inline void add_char_to_float_mult(float *arr, unsigned char* arr8,
    long nsamps, double scale)
{
#pragma omp parallel for
  for (int i=0; i<nsamps; i++)
    arr[i] = arr[i] + ((float) arr8[i])*scale;
}

static inline void float_to_char_div(float *arr, unsigned char* arr8,
    long nsamps, double div)
{
#pragma omp parallel for
  for (int i=0; i<nsamps; i++)
    arr8[i] = (unsigned char) round(arr[i]/div-0.01);
}

// 32bit
static inline void add_float_to_float_mult(float *arr, float* arr_in,
    long nsamps, double scale)
{
#pragma omp parallel for
  for (int i=0; i<nsamps; i++)
    arr[i] = arr[i] + arr_in[i]*scale;
}

static inline void float_to_float_div(float *arr, float* arr_in,
    long nsamps, double div)
{
#pragma omp parallel for
  for (int i=0; i<nsamps; i++)
    arr_in[i] = (float) round(arr[i]/div-0.01);
}



int main(int argc, char* argv[])
{
  char foutname[256];
  double weights[MAXFILES] = { [0 ... MAXFILES-1] = 1. };

  int iweight = 0;
  int arg = 0;
  int numthreads = 1;
  while ((arg=getopt(argc,argv,"o:w:p:")) != -1)
  {
    switch (arg)
    {
      case 'o':
        if (optarg)
        {
          if (sscanf (optarg, "%s", foutname) != 1)
          {
            fprintf(stderr, "ERROR: could not parse output filfile from %s\n",
                optarg);
            return EXIT_FAILURE;
          }
        }
        else
        {
          fprintf(stderr, "ERROR: -o requires argument\n");
          usage();
          return EXIT_FAILURE;
        }
        break;

      case 'w':
        if (optarg)
        {
          if (sscanf (optarg, "%lf", &weights[iweight]) != 1)
          {
            fprintf(stderr, "ERROR: could not parse in basedir from %s\n",
                optarg);
            return EXIT_FAILURE;
          }
          iweight += 1;
        }
        break;
      case 'p':
        if (optarg)
        {
          if (sscanf (optarg, "%i", &numthreads) != 1)
          {
            fprintf(stderr, "ERROR: could not parse in numthreads from %s\n",
                optarg);
            return EXIT_FAILURE;
          }
        }
    }
  }

  omp_set_num_threads(numthreads);

  int nargs = argc - optind;
  int ninfiles = nargs;

  // Only 1 file given, just cp it
  if (nargs == 1)
  {
    char cmd[64];
    sprintf(cmd, "cp %s %s", argv[optind], foutname);
    fprintf(stderr, "cmd: %s\n", cmd);
    int ret = system(cmd);
    if (ret == -1)
    {
      fprintf(stderr, "Could not execute: %s", cmd);
      return EXIT_FAILURE;
    }
    else
      return EXIT_SUCCESS;
  }
  else if (nargs <= 0)
  {
    fprintf(stderr, "ERROR: at least 1 fil files should be passed\n");
    usage();
    exit(EXIT_FAILURE);
  }

  if ((iweight !=0) & (iweight != ninfiles))
  {
    fprintf(stderr, "ERROR: number of weights (%i) is different than the "
        "number of input files (%i)\n", iweight, ninfiles);
    exit(EXIT_FAILURE);
  }

  char infnames[MAXFILES][256];

  for (arg=0; arg<nargs; arg++)
  {
    if (sscanf (argv[optind + arg], "%s", infnames[arg]) != 1)
    {
      fprintf(stderr, "ERROR: could not parse in basedir from %s\n",
          optarg);
      return EXIT_FAILURE;
    }
  }
  

  fil_t ** infilbanks;
  infilbanks = malloc(ninfiles * sizeof *infilbanks);

  int err=0;
  for (int i=0; i<nargs; i++)
  {
    infilbanks[i] = create_fil(infnames[i], &err, 
        SIGPROC_CREATE_FIL_READ);
    if (err < 0)
    {
      fprintf(stderr,"ERROR: could not create_fil on"
          " fname: %s\n", infnames[i]);
      return EXIT_FAILURE;
    }
  }

  fil_t * outfilbank = create_fil(foutname, &err, 
      SIGPROC_CREATE_FIL_WRITE);
  if (err < 0)
  {
    fprintf(stderr, "ERROR: could not create_fil for: %s\n",
        foutname);
    return EXIT_FAILURE;
  }

  // basic sanity checks
  if (sanity_check_files(infilbanks, ninfiles) == EXIT_FAILURE)
    return EXIT_FAILURE;

  // Write output header
  memcpy(outfilbank->header, infilbanks[0]->header,
      sizeof * outfilbank->header);
  write_header(*outfilbank);

  int nsamps_per_gulp = NSEC/infilbanks[0]->header->tsamp;
  unsigned long total_samps = 1e9;
  for (int i=0; i<ninfiles; i++)
    total_samps = MIN(total_samps, infilbanks[i]->header->nsamples);
  int nblocks = total_samps/nsamps_per_gulp;

  int nbit = infilbanks[0]->header->nbits;

  float *block;
  block = malloc(nsamps_per_gulp * infilbanks[0]->header->nchans
      *sizeof *block);
  unsigned char *block8bit;
  block8bit = malloc(nsamps_per_gulp * infilbanks[0]->header->nchans
      *sizeof *block8bit);

  double final_scale=0;
  for (int i=0; i<ninfiles; i++)
  {
    fseek(infilbanks[i]->file, infilbanks[i]->header->totalbytes, SEEK_SET);
    final_scale += weights[i];
  }


  size_t nmem;
  for (int iblock=0; iblock<nblocks+1; iblock++)
  {
    //clock_t start,end;
    //double cpu_time_used;
    //start = clock();

    int nsamps = MIN(nsamps_per_gulp, 
        total_samps - nsamps_per_gulp*iblock);
    long nsamps_chan = nsamps * infilbanks[0]->header->nchans;

    memset(block, 0, nsamps_chan*sizeof *block);

    for (int ifile=0; ifile<ninfiles; ifile++)
    {
      nmem = fread(block8bit, 1, nsamps_chan, infilbanks[ifile]->file);
      if (nmem != nsamps_chan)
      {
        fprintf(stderr,"ERROR on reading filterbank files: "
            "expected %li nsamps_chan, got %li. Exiting...\n",
            nsamps_chan, nmem);
        return EXIT_FAILURE;
      }

      add_char_to_float_mult(block, block8bit, nsamps_chan, weights[ifile]);
    }

    float_to_char_div(block, block8bit, nsamps_chan, final_scale);

    nmem = fwrite(block8bit, 1, nsamps_chan, outfilbank->file);
    if (nmem != nsamps_chan)
    {
      fprintf(stderr,"ERROR on writting output filterbank file, exiting...\n");
      return EXIT_FAILURE;
    }
    //end = clock();
    //cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;
    //fprintf(stderr,"Time: %lf\n", cpu_time_used);
  }




  // Destroy outputfilbank
  if (destroy_fil(outfilbank) < 0)
  {
    fprintf(stderr,"ERROR: could not destroy outfilbank\n");
    exit(EXIT_FAILURE);
  }

  // Destroy input filbanks
  for (int i=0; i<nargs; i++)
  {
    if (destroy_fil(infilbanks[i]) < 0)
    {
      fprintf(stderr,"ERROR: could not destroy outfilbank\n");
      exit(EXIT_FAILURE);
    }
  }

  free(infilbanks);
  free(block);
  free(block8bit);

  return EXIT_SUCCESS;
}
