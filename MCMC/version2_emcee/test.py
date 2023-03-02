import multiprocessing
import logging
import random

def addition_func(logger):
    a = random.randint(1,100)
    b = random.randint(1,100)
    c = a + b
    logger.info("Adding two random numbers: {} + {} = {}".format(a,b,c))

def worker(num):
    # Set up a logger with a unique output filename
    logger = logging.getLogger("Process {}".format(num))
    logger.setLevel(logging.INFO)
    file_handler = logging.FileHandler("process_{}.log".format(multiprocessing.current_process().pid))
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Log some information
    logger.info("Starting process {}".format(num))
    addition_func(logger)
    logger.info("Ending process {}".format(num))

# Create a list of processes
processes = []
for num in range(2):
    p = multiprocessing.Process(target=worker, args=(num,))
    processes.append(p)
    p.start()

# Wait for all processes to finish
for p in processes:
    p.join()
