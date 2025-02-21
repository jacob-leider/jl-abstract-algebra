#!/usr/bin/python3

# 

class LFSR:
    def __init__(self, coeffs: list, seed: list | None=None):
        """
        An inefficient linear feedback shift register.

            - output bit: 0
        """
        self._coeffs = coeffs
        if seed is None:
            self._seed = [0] * len(self._coeffs)
        else:
            self._seed = seed
        self._state = seed

    def set_seed(self, seed=None) -> list:
        

    def step(self):
        """
        Move forward one state.
        """
        # Shift indices downard.
        self._coeffs[:-1] = self._coeffs[1:]
        f = 0
        for i in range(len(self._coeffs) - 1):
            f += self._coeffs[i] * self._state[i]
        self._coeffs[-1] = f


    def state(self):
        return self._state


    def coeffs(self):
        return self._coeffs
