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

        return x[min_energy_i], E[min_energy_i]

    def plot(self, x, E):
    
        plt.plot(x, E, label=self.label)

        plt.xlabel('Bond length (a.u.)')

        plt.ylabel('Total Energy (a.u.)')

        plt.legend()

        plt.title(self.title)

        plt.show()
