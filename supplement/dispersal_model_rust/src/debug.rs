impl std::fmt::Debug for crate::Patch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Patch")
            .field("r", &self.resources)
            .field("/", &self.max_resources)
            .finish()
    }
}

impl std::fmt::Debug for crate::Family {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let geo = crate::hexgrid::geo_coordinates(self.location);

        f.debug_struct("Family")
            .field("descendence", &self.descendence)
            .field("size", &self.effective_size)
            .field("seasons_till_next_child", &self.seasons_till_next_child)
            .field("location", &self.location)
            .field("lat", &(geo.lat * 180. / crate::PI))
            .field("lon", &(geo.lon * 180. / crate::PI))
            .field("culture", &format_args!("{:020b}", &self.culture))
            .field("stored_resources", &self.stored_resources)
            .finish()
    }
}
