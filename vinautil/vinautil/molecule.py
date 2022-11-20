from pathlib import Path, PurePath
from openbabel import pybel
from vinautil.vinaconfig import xyz_point

from vinautil.utils.typecheck import typeassert


class center(xyz_point):
    def __init__(self, x, y, z):
        super().__init__(x, y, z)

    def __repr__(self):
        return 'center: {},{},{}'.format(self.x, self.y, self.z)

class box_size(xyz_point):
    def __init__(self, x, y, z):
        super().__init__(x, y, z)

    def __repr__(self):
        return 'box_size: {},{},{}'.format(self.x, self.y, self.z)


@typeassert(file=Path, fmt=str)
class molecule:

    def __init__(self, file: Path, fmt:str=None):
        self.file = file
        self.fmt = fmt if bool(fmt) else file.suffix.replace('.','')
        # self.center, self.box_size = self.compute_box().values()

    def compute_box(self, extending = 10):
        coords_lst = self.get_coord_lst()
        xlst, ylst, zlst = [i[0] for i in coords_lst], [i[1] for i in coords_lst], [i[2] for i in coords_lst]
        ([minX, minY, minZ], [maxX, maxY, maxZ]) = ([min(xlst), min(ylst), min(zlst)], [max(xlst), max(ylst), max(zlst)])
        minX = minX - float(extending)
        minY = minY - float(extending)
        minZ = minZ - float(extending)
        maxX = maxX + float(extending)
        maxY = maxY + float(extending)
        maxZ = maxZ + float(extending)
        SizeX = maxX - minX
        SizeY = maxY - minY
        SizeZ = maxZ - minZ
        CenterX = (maxX + minX) / 2
        CenterY = (maxY + minY) / 2
        CenterZ = (maxZ + minZ) / 2
        c = {
                'x': float("{:.2f}".format(CenterX)),
                'y': float("{:.2f}".format(CenterY)),
                'z': float("{:.2f}".format(CenterZ)),
        }
        s = {
                'x': float("{:.2f}".format(SizeX)),
                'y': float("{:.2f}".format(SizeY)),
                'z': float("{:.2f}".format(SizeZ))
        }
        rd = {
            'center': center(**c),
            'box_size': box_size(**s)
        }
        return rd

    @property
    def center(self):
        return self.center

    @property
    def box_size(self):
        return self.box_size

    def get_coord_lst(self, addHs:bool = True):
        molH = pybel.readfile(format = self.fmt, filename=self.file.__str__())
        molH = next(molH)
        # molH.OBMol.DeleteHydrogens()
        if addHs: molH.OBMol.AddPolarHydrogens()
        return [list(atom.coords) for atom in molH]

    def get_coord_dict(self, addHs:bool = True):
        molH = pybel.readfile(format = self.fmt, filename=self.file.__str__())
        molH = next(molH)
        if addHs: molH.OBMol.AddPolarHydrogens()
        return {atom.idx: atom.coords for atom in molH}

if __name__ == '__main__':
    file = Path('./test_file/2xd1_C_CEF.mol2')
    mol = molecule(file)