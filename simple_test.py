#!/usr/bin/env python3
"""
Simple test to verify the base score component fix.
"""

import os
import sys
import numpy as np
from typing import List

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent


class SimpleTestComponent(BaseScoreComponent):
    """Simple test component to verify step parameter passing."""
    
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.received_step = None
    
    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        """Implementation that records the step parameter."""
        print(f"calculate_score called with step={step}")
        self.received_step = step
        scores = np.array([0.5] * len(molecules), dtype=np.float32)
        return ComponentSummary(total_score=scores, parameters=self.parameters)


def main():
    print("🔧 Simple test for base score component fix")
    print("=" * 50)
    
    # Create test component
    params = ComponentParameters(
        component_type="test",
        name="Test",
        weight=1.0
    )
    component = SimpleTestComponent(params)
    
    # Test molecules
    molecules = ["CCO", "CCC"]
    test_step = 99
    
    print(f"\n1. Testing calculate_score directly with step={test_step}")
    result1 = component.calculate_score(molecules, step=test_step)
    print(f"   Received step: {component.received_step}")
    assert component.received_step == test_step, f"Direct call failed: {component.received_step}"
    
    print(f"\n2. Testing calculate_score_for_step with step={test_step}")
    component.received_step = None  # Reset
    result2 = component.calculate_score_for_step(molecules, step=test_step)
    print(f"   Received step: {component.received_step}")
    
    if component.received_step == test_step:
        print("✅ SUCCESS: Base class fix is working!")
        print(f"   calculate_score_for_step properly passed step={test_step}")
    else:
        print("❌ FAILURE: Base class fix is not working!")
        print(f"   Expected step={test_step}, got step={component.received_step}")
        return False
    
    print(f"\n3. Testing with default step (-1)")
    component.received_step = None  # Reset
    result3 = component.calculate_score_for_step(molecules)  # No step parameter
    print(f"   Received step: {component.received_step}")
    assert component.received_step == -1, f"Default step failed: {component.received_step}"
    
    print("\n✅ All tests passed!")
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
