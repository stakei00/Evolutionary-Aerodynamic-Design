class Airfoil:

    def __init__(self) -> None:
        """
        Initializes airfoil object 
        """
        self.chord = None 
        self.thickness_to_chord = None 
        self.max_camber_to_chord = None
        self.max_camber_from_le = None 

    def run_XFOIL(self) -> None: 
        """
        Analysis of airfoil in XFOIL. Adds results to object attributes
        """