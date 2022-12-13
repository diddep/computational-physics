#import runpy
import alpha_plot
import derivative_energy_plot
import energy_plot
import histogram_plot
import plot_task1_xdist
import plots_task2
import time

def main():
    start_time = time.time()

    alpha_plot.main()
    alpha_time = time.time()
    print("Alpha done in: %s seconds" % (alpha_time - start_time))

    derivative_energy_plot.main()
    derivative_time = time.time()
    print("Derivative done in: %s seconds" % (derivative_time - alpha_time))
    
    energy_plot.main()
    energy_time = time.time()
    print("Energy done in: %s seconds" % (energy_time - derivative_time))
    
    histogram_plot.main()
    histogram_time = time.time()
    
    print("Histogram done in: %s seconds" % (histogram_time - energy_time))
    
    plot_task1_xdist.main()
    task1_time = time.time()
    print("Task1 done in: %s seconds" % (task1_time - histogram_time))

    plots_task2.main()
    task2_time = time.time()
    print("Task2 done in: %s seconds" % (task2_time - task1_time))
    print("Finished all plots in: %s seconds" % (task2_time - task1_time))

main()

