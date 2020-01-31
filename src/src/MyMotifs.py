from MySeq import MySeq


class MyMotifs:

    def __init__(self, seqs):
        self.seqs = seqs.upper()
        self.tam = len(seqs[0])
        self.alfabeto = seqs[0].alfabeto()
        self.pwm = pwm

    def calculaContagens(self):
        return None

    def criaPWM(self):
        return None

    def consenso(self):
        return None

    def probabSeq(self):
        return None

    def seqMaisProvavel(self, seq):
        return None