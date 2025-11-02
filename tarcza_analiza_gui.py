from numpy import set_printoptions
from gui import MainApplication

set_printoptions(precision=15, suppress=False)

if __name__ == "__main__":
    app = MainApplication()
    app.mainloop()