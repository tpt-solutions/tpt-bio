# HCS test fixtures

Metadata written by the reference OME-NGFF writer, `ome-zarr` 0.21.0 (BSD-3-Clause), used
as a test oracle only. Its functions `write_plate_metadata`, `write_well_metadata` and
`write_image` made a 2-row, 3-column plate with two acquisitions, two wells (`A/1`, `B/2`)
and two fields per well. Each field is a 3-channel (`c`, `y`, `x`) 8x8 image. Channel labels
were added to each field's `omero` block afterwards, because the writer needs an `Omero`
object. The pixel data is synthetic and is not copied here, only the JSON metadata.

- `0.4/`: `.zattrs` files, NGFF 0.4 layout.
- `0.5/`: `zarr.json` files, NGFF 0.5 layout (metadata under `attributes.ome`).
