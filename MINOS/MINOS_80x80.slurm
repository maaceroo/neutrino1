#!/bin/bash
#
#SBATCH --job-name=MINOS
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --partition=icn

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=alexis@nucleares.unam.mx,marioacero@mail.uniatlantico.edu.co

# setup system
source load-modules.sh

echo -e "\n## Slurm job started on $(date +'%d-%m-%Y at %T') #####################"

#-----------------------------------------------------------------------------
# Inicia trabajo

export JOBNUM=`echo $SLURM_JOB_ID`

mkdir ${JOBNUM}
mkdir ${JOBNUM}/data
mkdir ${JOBNUM}/files_plots

export FUDGEE=1.00
export FUDGEN=1.00
export FUDGEBE=1.00
export FUDGEBN=1.00
export FUDGEBWSE=1.00
export FUDGEBWSN=1.00

export ENSCALE=1.00

srun ./MINOS_osc_ana_80x80.sh

# Termina trabajo
#-----------------------------------------------------------------------------
echo -e "## Slurm job finished on $(date +'%d-%m-%Y at %T') ####################"

#-----------------------------------------------------------------------------
# ouptut to slurm log file

echo "## Job name:                         $SLURM_JOB_NAME"
echo "## Job ID:                           $SLURM_JOB_ID"

