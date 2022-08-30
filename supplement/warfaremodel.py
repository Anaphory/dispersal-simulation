"""A toy model for warfare between different groups."""
import random
import collections
from dataclasses import dataclass

@dataclass
class F():
    c: int
    n: int

    def __str__(self):
        return ["#", "0", "$", "Δ"][self.c] * self.n

    __repr__ = __str__

    def __lt__(self, other):
        return (self.c, self.n) < (other.c, other.n)

def fight(families):
    # They come in random order.
    random.shuffle(families)

    print('So the population trying to move into the patch looks like this. (Not in that order.)')
    for f in sorted(families):
        print(f)

    groups = []
    while families:
        f = families.pop()
        print(f'In comes I, {f}.', end=" ")
        for group in groups:
            join = True
            for fx in group:
                print(f"I meet {fx}.", end=" ")
                if fx.c != f.c:
                    case = random.randrange(4)
                    if case == 0:
                        print(f"It's a massacre.", end=" ")
                        f.n, fx.n = max(f.n -fx.n, 0), max(fx.n -f.n, 0)
                    elif case == 1:
                        print(f"They ambush me.", end=" ")
                        f.n -= min(f.n, fx.n)
                    elif case == 2:
                        print(f"I slaughter them.", end=" ")
                        fx.n -= min(fx.n, f.n)
                    else:
                        print(f"We give each other a wide berth.", end=" ")
                        pass
                    join = False
                if f.n == 0:
                    print("I die.")
                    break
            else:
                if join:
                    print("I join.")
                    group.append(f)
                    break
                continue
            break
        else:
            if f.n > 0:
                print("I found a new tribe.")
                groups.append([f])
        groups = [[fx for fx in group if fx.n>0]
                for group in groups
                if [fx for fx in group if fx.n>0]]
        print(groups)
    return groups


# Okay, so there's a bunch of families of different size, each with some culture or another.
families = [
    F(int(random.randrange(16) ** 0.5), random.randrange(2, 11))
    for _ in range(random.randrange(1, 21))
]

fight(families)

c = collections.Counter()

for i in range(200):
    even = [F(0, 2), F(0, 3), F(0, 4), F(0, 5), F(0, 6)] + [F(1, 10) for _ in range(2)]
    groups = fight(even)
    for g in groups:
        for f in g:
            c[f.c] += f.n
    print(c)
