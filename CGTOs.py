import basis

def getH1s():

    H1s = basis.CGTO((0, 0, 0))
    H1s.add([3.425250910, 0.154328970])
    H1s.add([0.623913730, 0.535328140])
    H1s.add([0.168855400, 0.444634540])

    return H1s

def getHe1s():

    He1s = basis.CGTO((0, 0, 0))
    He1s.add([6.36242139*1.24**2, 0.154328970])
    He1s.add([1.15892300*1.24**2, 0.535328140])
    He1s.add([0.31364979*1.24**2, 0.444634540])

    return He1s

def getF1s():

    F1s = basis.CGTO((0, 0, 0))
    F1s.add([166.6791300, 0.15432897])     
    F1s.add([30.3608120, 0.53532814])     
    F1s.add([8.2168207, 0.44463454])

    return F1s

def getF2s():

    F2s = basis.CGTO((0, 0, 0))
    F2s.add([6.4648032, -0.09996723])
    F2s.add([1.5022812, 0.39951283])
    F2s.add([0.4885885, 0.70011547])

    return F2s

def getF2px():

    F2px = basis.CGTO((1, 0, 0))
    F2px.add([6.4648032, 0.15591627])
    F2px.add([1.5022812, 0.60768372])
    F2px.add([0.4885885, 0.39195739])

    return F2px

def getF2py():

    F2py = basis.CGTO((0, 1, 0))
    F2py.add([6.4648032, 0.15591627])
    F2py.add([1.5022812, 0.60768372])
    F2py.add([0.4885885, 0.39195739])

    return F2py

def getF2pz():

    F2pz = basis.CGTO((0, 0, 1))
    F2pz.add([6.4648032, 0.15591627])
    F2pz.add([1.5022812, 0.60768372])
    F2pz.add([0.4885885, 0.39195739])

    return F2pz
