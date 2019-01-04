import numpy as np

class CGTO():

    def __init__(self, powers):

        self.powers = powers
        self.infoList = []

    def add(self, info):

        self.infoList.append(info)

class basis():

    def __init__(self, CGTO):

        self.CGTO = CGTO

    def setCoords(self, N_coords):
        
        self.N_coords = N_coords
