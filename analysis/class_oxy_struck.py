__author__ = 'rasmus'


class OxyStruck:
    def __init__(self, oxygen_id, frame):
        self.in_bulk = False
        self.has_hbond = False
        self.frame = frame
        self.Hbonds = []
        self.oxygenID = oxygen_id

    def add_hbond(self, donor, acceptor, hydrogen):
        self.Hbonds.append(Hbond(donor, acceptor, hydrogen))


class Hbond:
    def __init__(self, donor, acceptor, hydrogen):
        self.donor = donor
        self.acceptor = acceptor
        self.hydrogen = hydrogen
        self.broken = False

    def __repr__(self):
        return str('d' + str(self.donor) + 'a' + str(self.acceptor) + 'h' + str(self.hydrogen))

    def __eq__(self, other):
        return self.donor == other.donor and self.acceptor == other.acceptor and self.hydrogen == other.hydrogen


class HbondStrcture:
    def __init__(self, hbond, start):
        self.framelist = []
        self.framelist.append(start)
        self.hbond = hbond

    def __repr__(self):
        return str(self.hbond)
