#!/usr/bin/env python3
# Name: Ryan Wanless (rwanless)
# Group Members: Boxuan Ma, Jesse Smith, Shreya Sinha

'''
Program docstring goes here
'''

import math


class Triad:
    """
    Calculate angles and distances among a triad of points.

    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()

    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r

    """

    def __init__(self, p, q, r):
        """ Construct a Triad.

        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ).
        """
        self.p = p
        self.q = q
        self.r = r

    # private helper methods
    def d2(self, a, b):  # calculate squared distance of point a to b
        return float(sum((ia - ib) * (ia - ib) for ia, ib in zip(a, b)))

    def dot(self, a, b):  # dotProd of standard vectors a,b
        return float(sum(ia * ib for ia, ib in zip(a, b)))

    def ndot(self, a, b, c):  # dotProd of vec. a,c standardized to b
        return float(sum((ia - ib) * (ic - ib) for ia, ib, ic in zip(a, b, c)))

    # calculate lengths(distances) of segments PQ, PR and QR
    def dPQ(self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p, self.q))

    def dPR(self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p, self.r))

    def dQR(self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q, self.r))

    def angleP(self):
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q, self.p, self.r) /
                         math.sqrt(self.d2(self.q, self.p) * self.d2(self.r, self.p)))

    def angleQ(self):
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p, self.q, self.r) /
                         math.sqrt(self.d2(self.p, self.q) * self.d2(self.r, self.q)))

    def angleR(self):
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p, self.r, self.q) /
                         math.sqrt(self.d2(self.p, self.r) * self.d2(self.q, self.r)))


def main():
    #gets input from uses
    coordinateInput = input("please put cordinates here")
    #replaces ( and ) with a , and splits via ,
    CORD = coordinateInput.replace('(', ',').replace(')', ',').split(',')
    #prints out corresponding stings from list and converts them into float
    p = (float(CORD[1]), float(CORD[2]), float(CORD[3]))
    q = (float(CORD[5]), float(CORD[6]), float(CORD[7]))
    r = (float(CORD[9]), float(CORD[10]), float(CORD[11]))
    #puts values of p q and r into class triad and renames it to mytriad
    myTriad = Triad(p, q, r)
    #runs triad with numbers from input to method dPQ
    pqlen = myTriad.dPQ()
    prlen = myTriad.dPR()
    #runs triad with numbers from input to meth dQR
    qrlen = myTriad.dQR()
    #get angle of Q and converts it to degres
    angleQ = myTriad.angleQ() * ((180 / math.pi))
    #prints out corresponding bond values
    print('N-C bond length:{0:0.2f}'.format(pqlen))
    print('N-Ca bond length:{0:0.2f}'.format(qrlen))
    print('N-Ca bond angle:{0:0.1f}'.format(angleQ))

    pass


main()