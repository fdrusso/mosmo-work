import collections
import copy

import chemparse
from ipydatagrid import DataGrid
from mosmo.model import DS, DbXref, Molecule
import pandas as pd


class MolEditTask:
    FOCUS_DBS = [DS.KEGG, DS.METACYC, DS.CHEBI, DS.BIGG]

    def __init__(self, ref_mol, kb_mol):
        self.ref_mols = [ref_mol]
        self.kb_mol = kb_mol
        self.dirty = False
        self.defunct_kb_id = None

    def data(self):
        data_mol = self.kb_mol
        data = {
            'ref_id': ', '.join(ref_mol.id for ref_mol in self.ref_mols),
            'kb_id': self.kb_mol.id or '',
            'name': data_mol.name,
            'shorthand': data_mol.shorthand or '',
        }

        xrefs = {xref.db: xref.id for xref in (data_mol.xrefs or [])}
        for db in MolEditTask.FOCUS_DBS:
            if db in xrefs:
                data[db.id] = xrefs.pop(db)
            else:
                data[db.id] = ''

        others = sorted(xrefs.items(), key=lambda kv: kv[0].id)
        data['others'] = ', '.join(f'{db.id}:{id}' for db, id in others)

        return data

    def update_kb_id(self, id):
        if not id:
            return
        if id != self.kb_mol.id:
            if self.kb_mol.id is not None:
                self.defunct_kb_id = self.kb_mol.id
            self.kb_mol.id = id
            self.dirty = True

    def update_name(self, name):
        if not name:
            return
        if name != self.kb_mol.name:
            # Preserve the old name as an aka, and remove the new one if it was an aka.
            if self.kb_mol.aka is not None:
                self.kb_mol.aka.insert(0, self.kb_mol.name)
                try:
                    self.kb_mol.aka.remove(name)
                except ValueError:
                    pass
            else:
                self.kb_mol.aka = [self.kb_mol.name]

            self.kb_mol.name = name
            self.dirty = True

    def update_shorthand(self, shorthand):
        shorthand = shorthand or None
        if shorthand != self.kb_mol.shorthand:
            self.kb_mol.shorthand = shorthand
            self.dirty = True

    def update_xref(self, db, id):
        if self.kb_mol.xrefs is None:
            self.kb_mol.xrefs = set()

        existing = None
        for xref in self.kb_mol.xrefs:
            if xref.db == db:
                existing = xref
                break
        if existing is not None and existing.id == id:
            return

        if existing is not None:
            self.kb_mol.xrefs.remove(existing)
            self.dirty = True
        if id:
            self.kb_mol.xrefs.add(DbXref(db=db, id=id))
            self.dirty = True

    def update_other_xrefs(self, xrefs_str):
        xrefs = set()
        if xrefs_str:
            for xref_str in xrefs_str.split(','):
                db, id = xref_str.strip().split(':')
                xrefs.add(DbXref(db=DS.get(db), id=id))
        for xref in xrefs:
            if xref not in self.kb_mol.xrefs:
                self.dirty = True
        for xref in self.kb_mol.xrefs:
            if xref.db in MolEditTask.FOCUS_DBS:
                xrefs.add(xref)
            elif xref not in xrefs:
                self.dirty = True
        self.kb_mol.xrefs = xrefs


class MolEditGrid:
    def __init__(self, tasks):
        self.tasks = tasks
        df = pd.DataFrame([task.data() for task in self.tasks]).set_index('ref_id')
        self.grid = DataGrid(df, editable=True)
        self.grid.on_cell_change(self.on_cell_change)

    def on_cell_change(self, cell):
        # The most robust way to do this is to apply the entire row's data every time
        task = self.tasks[cell['row']]
        row = self.grid.data.iloc[cell['row']]

        task.update_kb_id(row['kb_id'])
        task.update_name(row['name'])
        task.update_shorthand(row['shorthand'])
        for db in task.FOCUS_DBS:
            task.update_xref(db, row[db.id])
        task.update_other_xrefs(row['others'])


class RxnEditTask:
    FOCUS_DBS = [DS.EC, DS.KEGG, DS.METACYC, DS.RHEA, DS.BIGG]

    def __init__(self, ref_rxn, kb_rxn):
        self.ref_rxn = ref_rxn
        self.kb_rxn = kb_rxn
        self.dirty = False
        self.defunct_kb_id = None

    def data(self):
        data = {
            'ref_id': self.ref_rxn.id,
            'kb_id': self.kb_rxn.id or '',
            'name': self.kb_rxn.name,
            'shorthand': self.kb_rxn.shorthand or '',
            'catalyst': self.kb_rxn.catalyst.id if self.kb_rxn.catalyst else '',
            'rev': self.kb_rxn.reversible,
            'equation': self.kb_rxn.equation,
        }

        xrefs = {xref.db: xref.id for xref in (self.kb_rxn.xrefs or [])}
        for db in RxnEditTask.FOCUS_DBS:
            if db in xrefs:
                data[db.id] = xrefs.pop(db)
            else:
                data[db.id] = ''

        others = sorted(xrefs.items(), key=lambda kv: kv[0].id)
        data['others'] = ', '.join(f'{db.id}:{id}' for db, id in others)

        return data

    def update_kb_id(self, id):
        if not id:
            return
        if id != self.kb_rxn.id:
            if self.kb_rxn.id is not None:
                self.defunct_kb_id = self.kb_rxn.id
            self.kb_rxn.id = id
            self.dirty = True

    def update_name(self, name):
        if not name:
            return
        if name != self.kb_rxn.name:
            # Preserve the old name as an aka, and remove the new one if it was an aka.
            if self.kb_rxn.aka is not None:
                self.kb_rxn.aka.insert(0, self.kb_rxn.name)
                try:
                    self.kb_rxn.aka.remove(name)
                except ValueError:
                    pass
            else:
                self.kb_rxn.aka = [self.kb_rxn.name]

            self.kb_rxn.name = name
            self.dirty = True

    def update_shorthand(self, shorthand):
        shorthand = shorthand or None
        if shorthand != self.kb_rxn.shorthand:
            self.kb_rxn.shorthand = shorthand
            self.dirty = True

    def update_catalyst(self, catalyst_id):
        if self.kb_rxn.catalyst:
            if self.kb_rxn.catalyst.id == id:
                return
            else:
                self.kb_rxn.catalyst = None
                self.dirty = True

        if catalyst_id:
            self.kb_rxn.catalyst = Molecule(id=catalyst_id)
            self.dirty = True

    def update_reversible(self, reversible):
        if reversible != self.kb_rxn.reversible:
            self.kb_rxn.reversible = reversible
            self.dirty = True

    def update_xref(self, db, id):
        if self.kb_rxn.xrefs is None:
            self.kb_rxn.xrefs = set()

        existing = None
        for xref in self.kb_rxn.xrefs:
            if xref.db == db:
                existing = xref
                break
        if existing is not None and existing.id == id:
            return

        if existing is not None:
            self.kb_rxn.xrefs.remove(existing)
            self.dirty = True
        if id:
            self.kb_rxn.xrefs.add(DbXref(db=db, id=id))
            self.dirty = True

    def update_other_xrefs(self, xrefs_str):
        xrefs = set()
        if xrefs_str:
            for xref_str in xrefs_str.split(','):
                db, id = xref_str.strip().split(':')
                xrefs.add(DbXref(db=DS.get(db), id=id))
        for xref in self.kb_rxn.xrefs:
            if xref.db in RxnEditTask.FOCUS_DBS:
                xrefs.add(xref)

        if xrefs != self.kb_rxn.xrefs:
            self.dirty = True
            self.kb_rxn.xrefs = xrefs


class RxnEditGrid:
    def __init__(self, tasks):
        self.tasks = tasks
        df = pd.DataFrame([task.data() for task in self.tasks]).set_index('ref_id')
        self.grid = DataGrid(df, editable=True)
        self.grid.on_cell_change(self.on_cell_change)

    def on_cell_change(self, cell):
        # The most robust way to do this is to apply the entire row's data every time
        task = self.tasks[cell['row']]
        row = self.grid.data.iloc[cell['row']]

        task.update_kb_id(row['kb_id'])
        task.update_name(row['name'])
        task.update_shorthand(row['shorthand'])
        task.update_catalyst(row['catalyst'])
        task.update_reversible(row['rev'])
        for db in task.FOCUS_DBS:
            task.update_xref(db, row[db.id])
        task.update_other_xrefs(row['others'])


def collect_reference_objects(ec_file, KB, ref_db):
    # Collect reference reactions cross-referenced to EC numbers
    ref_rxns = []
    with open(ec_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            ec_ref = DbXref(db=DS.EC, id=line)
            rxns = KB.xref(ref_db, ec_ref)
            if rxns:
                for rxn in rxns:
                    if rxn not in ref_rxns:
                        ref_rxns.append(rxn)
            else:
                print(f'{ec_ref} has no corresponding reactions in {ref_db.id}')

    # Collect reference molecules used in these reactions
    ref_mols = {ref_mol
                for rxn in ref_rxns
                for ref_mol in rxn.stoichiometry}

    return ref_rxns, ref_mols


def build_mol_edit_tasks(ref_mols, KB, kb_mol_db):
    def build_kb_mol(ref_mol):
        # Constructs a new kb molecule based on the ref, and cross-referenced to it.
        kb_mol = copy.deepcopy(ref_mol)
        kb_mol.id = None  # Must be assigned to save it.
        kb_mol.db = None  # Assigned by KB on save
        if kb_mol.xrefs:
            kb_mol.xrefs.add(ref_mol.ref)
        else:
            kb_mol.xrefs = {ref_mol.ref}
        return kb_mol

    # Keep track of tasks by kb_mol to be edited.
    # Because of variant forms, multiple ref mols may map to one canonical kb mol.
    tasks = []
    covered = {}
    for ref_mol in ref_mols:
        # Use existing kb molecules where possible.
        kb_mols = KB.xref(kb_mol_db, ref_mol.ref)
        if len(kb_mols) > 1:
            print(f'{ref_mol} maps to multiple KB entries.')
        if kb_mols:
            kb_mol = kb_mols[0]
            if kb_mol.canonical_form is not None:
                kb_mol = KB.get(kb_mol_db, kb_mol.canonical_form.parent_id)
        else:
            kb_mol = build_kb_mol(ref_mol)

        if kb_mol in covered:
            covered[kb_mol].ref_mols.append(ref_mol)
        else:
            task = MolEditTask(ref_mol, kb_mol)
            tasks.append(task)
            if kb_mol.id is not None:
                covered[kb_mol] = task

    # The original order of molecules is not meaningful. Sort them for easier editing.
    return sorted(tasks, key=lambda task: task.kb_mol.name)


def build_rxn_edit_tasks(ref_rxns, KB, kb_rxn_db, kb_mol_db):
    def build_kb_rxn(ref_rxn):
        # Constructs a new kb reaction based on the ref, and cross-referenced to it.
        # Reactions using molecules missing from the kb are skipped.
        stoichiometry = {}
        for ref_mol, count in ref_rxn.stoichiometry.items():
            kb_mols = KB.xref(kb_mol_db, ref_mol.ref)
            if len(kb_mols) == 1:
                kb_mol = kb_mols[0]
                if kb_mol.canonical_form is not None:
                    kb_mol = KB.get(kb_mol_db, kb_mol.canonical_form.parent_id)

                stoichiometry[kb_mol] = count
            else:
                # Abort
                return None

        kb_rxn = copy.deepcopy(ref_rxn)
        kb_rxn.id = None  # Must be assigned to save it.
        kb_rxn.db = None  # Assigned by KB on save
        kb_rxn.stoichiometry = stoichiometry
        if kb_rxn.xrefs is not None:
            kb_rxn.xrefs.add(ref_rxn.ref)
        else:
            kb_rxn.xrefs = {ref_rxn.ref}

        return kb_rxn

    tasks = []
    skipped = []
    for ref_rxn in ref_rxns:
        # Use existing kb reactions where possible.
        kb_rxns = KB.xref(kb_rxn_db, ref_rxn.ref)
        if len(kb_rxns) > 1:
            print(f'{ref_rxn} maps to multiple KB entries.')
        if kb_rxns:
            kb_rxn = kb_rxns[0]
        else:
            kb_rxn = build_kb_rxn(ref_rxn)

        if kb_rxn:
            tasks.append(RxnEditTask(ref_rxn, kb_rxn))
        else:
            skipped.append(ref_rxn)

    # Preserve the order of the reactions from the original input
    # return sorted(tasks, key=lambda task: task.kb_rxn.name), skipped
    return tasks, skipped

def mass_imbalance(reaction):
    net = collections.defaultdict(float)
    for mol, stoich in reaction.stoichiometry.items():
        for component, count in chemparse.parse_formula(mol.formula).items():
            net[component] += stoich * count
    return {component: count for component, count in net.items() if count != 0}

def charge_imbalance(reaction):
    net = 0
    for mol, stoich in reaction.stoichiometry.items():
        net += stoich * mol.charge
    return net

def check_reaction(reaction, ref_reaction=None):
    imbalance = mass_imbalance(reaction)
    charge = charge_imbalance(reaction)
    if imbalance or charge:
        print(f'{reaction}: mass imbalance = {imbalance}, charge imbalance = {charge}')
        print()
        print('reaction:')
        for mol, count in reaction.stoichiometry.items():
            print(f'  {count:>2d} {mol.formula:20} {mol}')
        if ref_reaction:
            print()
            print('reference:')
            for mol, count in ref_reaction.stoichiometry.items():
                print(f'  {count:>2d} {mol.formula:20} {mol}')
