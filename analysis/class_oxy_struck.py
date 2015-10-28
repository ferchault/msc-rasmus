__author__ = 'rasmus'


class OxyStruck:
    def __init__(self, oxygen_id, frame):
        self.in_bulk = False
        self.has_hbond = False
        self.frame = frame
        self.Hbonds = []
        self.oxygenID = oxygen_id

    def add_hbond(self, donor, acceptor, hydrogen):
        self.Hbonds.append(HbondStruck(donor, acceptor, hydrogen))


class HbondStruck:
    def __init__(self, donor, acceptor, hydrogen):
        self.donor = donor
        self.acceptor = acceptor
        self.hydrogen = hydrogen