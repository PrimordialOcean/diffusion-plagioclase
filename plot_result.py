import matplotlib.pyplot as plt
import pandas as pd

def main():
    times = [ 1, 100, 365, 365*20, 365*100, 100000]
    df = pd.read_csv('tmp.csv')
    
    df_ini = pd.read_csv("initial_value.csv")
    xan_measured = (df_ini["distance"], df_ini["XAn"])
    
    df_m = pd.read_csv('measured_value.csv')
    measured_data = (df_m["Distance(um)"],df_m["MgO"])
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    fig = plt.figure()
    plt.rcParams["font.size"] = 14
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ax.set_title("$T=1000 ^\circ$C (van Orman, 2014)")
    #ax.set_title(str(time_days)+" days in " + str(tempc) + "$^\circ$C")
    best_fit = 100000
    #ax.plot(*measured_data, "o", color="w", mec="k", label="Measured")
    for t in times:
        x = df.iloc[:,0]
        y = df.iloc[:,t]
        plotdata = (x, y)
        if t == 1:
            ax.plot(*plotdata, "--", color="k", label="Initial")
        elif t == best_fit:
            ax.plot(*plotdata, "-", color="r", label=str(int(t/365))+"yrs")
        else:
            if t < 365:
                ax.plot(x, y, "-", color="grey", label=str(int(t))+"days")
            else:
                ax.plot(x, y, "-", color="grey", label=str(int(t/365))+"yrs")
    ax2.plot(*xan_measured, "o", color="b", label="XAn")
    ax2.set_ylabel("XAn", fontsize=16)
    ax2.set_ylim(0.55, 2)
    ax.set_ylim(0, 0.7)
    ax.set_xlabel("Distance from rim (\u03bcm)", fontsize=16)
    ax.set_ylabel("MgO (wt.%)", fontsize=16)
    fig.legend(loc=1, fancybox=False, framealpha=1, edgecolor="k", fontsize=10)
    fig.savefig('img.jpg', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main()