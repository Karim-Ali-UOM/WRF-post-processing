from src.dictReader import readDict
from os import getcwd
from src.utilities import inputFileName, message

def main():
    message("Started WRF post-processing", 0)
    globe, ordersList = readDict(getcwd() + "/" + inputFileName())
    for order in ordersList: order.execute(globe)

main()