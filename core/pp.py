import matplotlib.pyplot as plt

class pp():

    def __init__(self, title, label):

        self.title = title
        self.label = label

    def min_energy(self, x, E):

        for i in range(0, len(x)):
            if i == 0:
                min_energy_i = 0
            else:
                if E[i] < E[min_energy_i]:
                    min_energy_i = i

        self.min_energy_i = min_energy_i
        self.x_min = x[min_energy_i]
        self.E_min = E[min_energy_i]

        return x[min_energy_i], E[min_energy_i]

    def plot(self, x, E):
    
        plt.plot(x, E, label=self.label, marker = 'o', markevery = [self.min_energy_i])

        plt.title(self.title)

        plt.xlabel('Bond length (a.u.)')
        plt.ylabel('Total Energy (a.u.)')

        x_min = self.x_min
        E_min = self.E_min

        text = '('+str(round(x_min, 3))+', '+str(round(E_min, 3))+')'
        plt.annotate(text, xy = (x_min, E_min), xytext = (x_min+1, E_min+1), arrowprops = dict(arrowstyle = '->'))

        plt.legend()

        plt.grid()
        plt.show()
