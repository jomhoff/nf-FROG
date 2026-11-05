import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HELPER = Path(__file__).resolve().parents[1] / 'bin/compare_references.py'

class ConvergenceTests(unittest.TestCase):
    def compare(self, previous, current, new, fraction=None, round_number=3, policy=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            task = root / 'work/aa/task'
            task.mkdir(parents=True)
            (task / '.command.sh').touch()
            for name, sequence in [('previous', previous), ('current', current), ('new', new)]:
                (task / (name + '.fa')).write_text('>chr1\n' + sequence + '\n')
            command = [sys.executable, str(HELPER), '--previous', 'previous.fa', '--current', 'current.fa',
                       '--new', 'new.fa', '--round', str(round_number), '--metrics', 'metrics.tsv',
                       '--status', 'status.txt']
            if fraction is not None:
                if policy:
                    (root / 'work/.frog-convergence-policy.json').write_text(json.dumps({
                        'maximum_change_fraction': fraction, 'consecutive_rounds': 2}))
                else:
                    command += ['--maximum-change-fraction', fraction]
            subprocess.run(command, cwd=task, check=True, capture_output=True)
            details = json.loads((task / 'metrics.details.json').read_text())
            self.assertEqual(len((task / 'metrics.tsv').read_text().splitlines()[1].split('\t')), 5)
            return (task / 'status.txt').read_text().strip(), details

    def test_strict_default_and_exact_flag(self):
        status, details = self.compare('A'*20000, 'C'+'A'*19999, 'C'+'A'*19999)
        self.assertEqual(status, 'converged')
        self.assertTrue(details['exact_zero_change_convergence'])

    def test_two_small_rounds_converge(self):
        status, details = self.compare('A'*20000, 'C'+'A'*19999, 'CC'+'A'*19998, '0.0001')
        self.assertEqual(status, 'converged')
        self.assertFalse(details['exact_zero_change_convergence'])
        self.assertEqual(details['previous_round_changed_bases'], 1)

    def test_large_previous_round_prevents_stopping(self):
        status, _ = self.compare('A'*20000, 'CC'+'A'*19998, 'CCC'+'A'*19997, '0.0001')
        self.assertEqual(status, 'active')

    def test_threshold_is_strict(self):
        status, _ = self.compare('A'*20000, 'C'+'A'*19999, 'CCC'+'A'*19997, '0.0001')
        self.assertEqual(status, 'active')

    def test_work_policy_and_oscillation_accounting(self):
        status, details = self.compare('A'*20000, 'C'+'A'*19999, 'A'*20000, '0.0001', policy=True)
        self.assertEqual(status, 'converged')
        self.assertEqual(details['oscillating_bases'], 1)
        self.assertTrue(details['policy_file'])
        strict_status, _ = self.compare('A'*20000, 'C'+'A'*19999, 'A'*20000)
        self.assertEqual(strict_status, 'active')

    def test_first_round_cannot_supply_two_rounds(self):
        status, _ = self.compare('A'*20000, 'A'*20000, 'A'*20000, '0.0001', round_number=1)
        self.assertEqual(status, 'active')

if __name__ == '__main__':
    unittest.main()
