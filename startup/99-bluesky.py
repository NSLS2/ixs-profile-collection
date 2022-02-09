'''
from bluesky.callbacks.olog import OlogCallback
from bluesky.callbacks import LivePlot


from suitcase.spec import DocumentToSpec

spec_cb = DocumentToSpec('/tmp/spec1.spec')
RE.subscribe(spec_cb)
'''


from suitcase.utils import MultiFileManager
from suitcase.specfile import Serializer
from event_model import RunRouter
import event_model
from pathlib import Path
from datetime import datetime


class MultiFileManagerHack(MultiFileManager):

    def open(self, label, postfix, mode, **kwargs):
        mode = 'a' if mode=='x' else mode
        f = super().open(label, postfix, mode, **kwargs)
        return f

    def reserve_name(self, label, postfix):
        name = (self._directory / Path(postfix)).expanduser().resolve()
        if name in self._reserved_names:
            return name
        return super().reserve_name(label, postfix)


class SerializerHack(Serializer):

    def event_page(self, doc):
        for event in event_model.unpack_event_page(doc):
            self.event(event)
        return doc

    def event(self, doc):
        doc = super().event(doc)
        self._file.flush()
        return doc


year, month, day = map(str, [datetime.now().year,
                             datetime.now().month,
                             datetime.now().day])

prefix = f'ixs_spec_{year}_{month}_{day}'
directory = '/home/xf10id/specfiles/'


def spec_factory(name, doc):
    spec_cb = SerializerHack(MultiFileManagerHack(directory, allowed_modes=('x','a')),
                     file_prefix=prefix)
    spec_cb(name, doc)
    return [spec_cb], []

run_router = RunRouter([spec_factory])
RE.subscribe(run_router)
