import os
import sys
from pathlib import Path

system = sys.argv[1]

directory=os.getcwd()
directory=directory+'/MCMC_'+system

combined_filename='AcceptedParameters.txt'
if os.path.exists(combined_filename):os.remove(combined_filename)
Path(combined_filename).touch()

for file in os.listdir(directory):
	filename=os.fsdecode(file)
	if filename.startswith('accepted'):
		accepted_filename=directory+'/'+filename
		with open(combined_filename,'a') as f1:
			with open(accepted_filename,'r') as f2:
				for i,lines in enumerate(f2):
					if i==0:
						continue
					if i>0:

						x=lines.split()
						current_state=x[1:-1]
						IterationNumber=float(x[0])
						if IterationNumber==1:
							Step=1
							previous_state=current_state
							f1.write(lines)
							continue
						else:
							Step=Step+1
							if IterationNumber!=Step:
								while True:
									if Step!=IterationNumber:
										state=[str(Step)]+previous_state
										Parameters='\t'.join(state)
										f1.write(Parameters+'\n')
										Step=Step+1
									else:
										state=[str(Step)]+current_state
										Parameters='\t'.join(state)
										f1.write(Parameters+'\n')
										previous_state=current_state
										break
							else:
								state=[str(Step)] + current_state
								Parameters='\t'.join(state)
								f1.write(Parameters+'\n')
								previous_state=current_state

