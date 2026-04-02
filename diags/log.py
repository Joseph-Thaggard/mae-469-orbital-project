import os

def log_start():
    if os.path.exists("output_log.txt"):
        os.remove("output_log.txt")
    with open("output_log.txt", "w") as log_file:
        log_file.write("Log started\n")

def log(string):
    with open("output_log.txt", "a") as log_file:
        log_file.write(string + "\n")