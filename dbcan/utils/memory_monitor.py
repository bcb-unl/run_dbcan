"""Memory monitoring utilities for dbcan."""
import logging
import psutil
import os
import time
from typing import Optional, Dict, List
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class MemoryMonitor:
    """Monitor memory usage and provide warnings."""
    
    def __init__(self, max_memory_usage: float = 0.8, check_interval: float = 1.0):
        """
        Initialize memory monitor.
        
        Args:
            max_memory_usage: Maximum memory usage ratio (0.0-1.0)
            check_interval: Interval between checks in seconds
        """
        self.max_memory_usage = max_memory_usage
        self.check_interval = check_interval
        self.process = psutil.Process(os.getpid())
        self._initial_memory = self.get_memory_info()
        self._memory_history: List[Dict] = []
        self._peak_memory: Optional[Dict] = None
        self._start_time: Optional[float] = None
    
    def get_memory_info(self) -> Dict[str, float]:
        """Get current memory information."""
        mem_info = self.process.memory_info()
        vm = psutil.virtual_memory()
        return {
            'rss_mb': mem_info.rss / (1024 * 1024),  # Resident Set Size in MB
            'vms_mb': mem_info.vms / (1024 * 1024),  # Virtual Memory Size in MB
            'available_mb': vm.available / (1024 * 1024),
            'total_mb': vm.total / (1024 * 1024),
            'percent': vm.percent,
            'used_mb': vm.used / (1024 * 1024),
        }
    
    def get_memory_usage_ratio(self) -> float:
        """Get current memory usage ratio (0.0-1.0)."""
        vm = psutil.virtual_memory()
        return vm.percent / 100.0
    
    def is_memory_safe(self) -> bool:
        """Check if memory usage is below threshold."""
        return self.get_memory_usage_ratio() < self.max_memory_usage
    
    def get_available_memory_mb(self) -> float:
        """Get available memory in MB."""
        vm = psutil.virtual_memory()
        return vm.available / (1024 * 1024)
    
    def log_memory_status(self, context: str = ""):
        """Log current memory status."""
        mem_info = self.get_memory_info()
        context_str = f" [{context}]" if context else ""
        logger.info(
            f"Memory status{context_str}: "
            f"RSS={mem_info['rss_mb']:.1f}MB, "
            f"Available={mem_info['available_mb']:.1f}MB, "
            f"Usage={mem_info['percent']:.1f}%"
        )
    
    def check_and_warn(self, context: str = "") -> bool:
        """
        Check memory usage and warn if approaching limit.
        
        Returns:
            True if memory is safe, False if approaching limit
        """
        usage_ratio = self.get_memory_usage_ratio()
        if usage_ratio >= self.max_memory_usage:
            mem_info = self.get_memory_info()
            context_str = f" in {context}" if context else ""
            logger.warning(
                f"Memory usage high{context_str}: {mem_info['percent']:.1f}% "
                f"(threshold: {self.max_memory_usage * 100:.1f}%). "
                f"Available: {mem_info['available_mb']:.1f}MB"
            )
            return False
        return True
    
    def estimate_batch_size(self, sequence_size_mb: float, safety_factor: float = 0.5) -> int:
        """
        Estimate safe batch size based on available memory.
        
        Args:
            sequence_size_mb: Size of one sequence in MB (estimated)
            safety_factor: Safety factor to leave memory headroom (0.0-1.0)
        
        Returns:
            Estimated number of sequences that can be processed in one batch
        """
        available_mb = self.get_available_memory_mb()
        # Reserve some memory for other operations
        usable_mb = available_mb * safety_factor
        
        if sequence_size_mb <= 0:
            # Default: assume 0.01 MB per sequence (average protein)
            sequence_size_mb = 0.01
        
        # Account for pyhmmer overhead (typically 3-5x sequence size)
        overhead_factor = 4.0
        batch_size = int((usable_mb / (sequence_size_mb * overhead_factor)))
        
        # Set reasonable bounds
        batch_size = max(100, min(batch_size, 50000))
        
        logger.debug(
            f"Estimated batch size: {batch_size} sequences "
            f"(available: {available_mb:.1f}MB, "
            f"seq_size: {sequence_size_mb:.3f}MB, "
            f"safety: {safety_factor})"
        )
        
        return batch_size
    
    def start_monitoring(self):
        """Start monitoring session."""
        self._start_time = time.time()
        self._memory_history = []
        self._peak_memory = None
        self._initial_memory = self.get_memory_info()
    
    def record_checkpoint(self, context: str = ""):
        """Record a memory checkpoint."""
        mem_info = self.get_memory_info()
        checkpoint = {
            'time': time.time(),
            'context': context,
            'rss_mb': mem_info['rss_mb'],
            'available_mb': mem_info['available_mb'],
            'percent': mem_info['percent']
        }
        self._memory_history.append(checkpoint)
        
        # Update peak memory
        if self._peak_memory is None or mem_info['rss_mb'] > self._peak_memory['rss_mb']:
            self._peak_memory = {
                'rss_mb': mem_info['rss_mb'],
                'available_mb': mem_info['available_mb'],
                'percent': mem_info['percent'],
                'context': context
            }
    
    def generate_report(self) -> Dict:
        """
        Generate a comprehensive memory usage report.
        
        Returns:
            Dictionary containing memory statistics
        """
        current_mem = self.get_memory_info()
        elapsed_time = time.time() - self._start_time if self._start_time else 0
        
        report = {
            'initial_memory_mb': self._initial_memory['rss_mb'],
            'current_memory_mb': current_mem['rss_mb'],
            'peak_memory_mb': self._peak_memory['rss_mb'] if self._peak_memory else current_mem['rss_mb'],
            'peak_memory_context': self._peak_memory['context'] if self._peak_memory else 'N/A',
            'memory_increase_mb': current_mem['rss_mb'] - self._initial_memory['rss_mb'],
            'available_memory_mb': current_mem['available_mb'],
            'total_memory_mb': current_mem['total_mb'],
            'memory_usage_percent': current_mem['percent'],
            'checkpoints': len(self._memory_history),
            'elapsed_time_seconds': elapsed_time
        }
        
        return report
    
    def log_report(self, context: str = ""):
        """Log a comprehensive memory usage report."""
        report = self.generate_report()
        context_str = f" [{context}]" if context else ""
        logger.info(
            f"Memory report{context_str}: "
            f"Initial={report['initial_memory_mb']:.1f}MB, "
            f"Current={report['current_memory_mb']:.1f}MB, "
            f"Peak={report['peak_memory_mb']:.1f}MB ({report['peak_memory_context']}), "
            f"Increase={report['memory_increase_mb']:.1f}MB, "
            f"Usage={report['memory_usage_percent']:.1f}%, "
            f"Checkpoints={report['checkpoints']}, "
            f"Elapsed={report['elapsed_time_seconds']:.1f}s"
        )
    
    @contextmanager
    def monitor(self, context: str = ""):
        """Context manager to monitor memory usage during an operation."""
        self.log_memory_status(f"Before {context}" if context else "Before operation")
        try:
            yield self
        finally:
            self.log_memory_status(f"After {context}" if context else "After operation")
            if not self.check_and_warn(context):
                logger.warning(f"Memory usage increased during {context}")


def get_memory_monitor(max_memory_usage: Optional[float] = None) -> MemoryMonitor:
    """
    Get a memory monitor instance.
    
    Args:
        max_memory_usage: Maximum memory usage ratio (default: 0.8)
    
    Returns:
        MemoryMonitor instance
    """
    if max_memory_usage is None:
        max_memory_usage = 0.8
    return MemoryMonitor(max_memory_usage=max_memory_usage)

