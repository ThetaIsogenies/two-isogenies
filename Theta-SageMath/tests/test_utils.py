from sage.all import EllipticCurve, random_prime, GF

def _random_prime(B):
    while True:
        p = random_prime(B)
        if p % 4 == 3:
            return p
        
def _random_field(B):
    p = _random_prime(B)
    return GF(p**2, name="i", modulus=[1,0,1])

def _random_supersingular_curve(F):
    E = EllipticCurve(F, [1,0])
    K = E.random_point()
    E_rand = E.isogeny(K, algorithm="factored").codomain()

    return E_rand.montgomery_model()

def random_supersingular_curve(B=1000):
    F = _random_field(B)
    return _random_supersingular_curve(F)

def random_supersingular_curves(B=1000):
    F = _random_field(B)
    E1 = _random_supersingular_curve(F)
    E2 = _random_supersingular_curve(F)

    return E1, E2