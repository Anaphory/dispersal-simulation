impl std::fmt::Debug for crate::Patch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Patch")
            //.field("r", &self.resources.iter().sum::<f32>())
            //.field("/", &self.resources.iter().sum::<f32>())
            .finish()
    }
}

impl std::fmt::Debug for crate::Family {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Family")
            .field("descendence", &self.descendence)
            .field("size", &self.effective_size)
            .field("location", &self.location)
            .field("culture", &format_args!("{:020b}", &self.culture))
            .field("stored_resources", &self.stored_resources)
            .field("last_harvest", &self.last_harvest)
            .finish()
    }
}

impl std::fmt::Debug for crate::Culture {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        self.in_memory.fmt(f)
    }
}
