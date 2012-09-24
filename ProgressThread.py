import threading

class ProgressThread(threading.Thread):
   def __init__(self, task_queue, progress_queue):
      self.task_queue = task_queue
      self.progress_queue = progress_queue

   def run(self):
      while (not self.task_queue.empty()):
         if (not self.progress_queue.empty()):
            progress = self.progress_queue.get()
